'''Functions and classes used to interface with .nib files as created by Jim
Kent's nibFrag and faToNib utilities.'''

import glob
import math
import os
import struct
import sys
import warnings
from cStringIO import StringIO
from collections import defaultdict as dd

from chipsequtil import reverse_complement, get_file_parts, BEDFile


# module fields
NOMASK,MASK,HARDMASK = range(3)


class NibException(Exception) : pass


def _nib_fd(nib) :
    '''Returns filename and file descriptor for nib, detecting whether it is a \
    path or fd appropriately'''

    # check to see if nib is a file or a string
    if isinstance(nib,file) :
        nib_fn = nib.name
        nib.seek(0)
        nib_f = nib
    elif isinstance(nib,str) :
        nib_fn = nib
        nib_f = open(nib,'rb')
    else :
        raise NibException('Incompatible .nib argument %s with type %s, needs to \
        be either <type \'file\'> or <type \'str\'>'%(str(nib),type(nib)))

    return nib_fn, nib_f


def get_nib(nib,start=0,end=-1,strand='+',mask=NOMASK,name=None,dbHeader=None,tbaHeader=None) :
    '''Return a (header,sequence) tuple representing this nibFrag record'''
    headers = get_nib_header_batch(nib,[(start,end,strand,name,dbHeader,tbaHeader),])
    seqs = get_nib_seq_batch(nib,[(start,end,strand)],mask)
    return headers[0], seqs[0]


def get_nib_batch(nib,queries,mask=NOMASK) :
    '''Batch interface for fetching fasta records.  Returns tuple of lists
    (headers,sequences)'''
    headers = get_nib_header_batch(nib,queries)
    seqs = get_nib_seq_batch(nib,[x[:3] for x in queries],mask=mask)
    return headers, seqs


def get_nib_seq(nib,start=0,end=-1,strand='+',mask=NOMASK) :
    '''Extract subsequence from .nib file like Jim Kent's nibFrag utility.
    Default behavior is to return the entire sequence.

    Extract the nucleotide substring defined by the closed interval [start,end]
    from the sequence found in *nib_fn*.  *mask* parameter has the following
    possible values:

    chipsequtil.nib.NOMASK -- masked positions are not indicated (default)
    chipsequtil.nib.MASK -- masked positions are capitalized, normal bases lower case
    chipsequtil.nib.NOMASK -- masked positions are replaced with Ns
    '''
    return get_nib_seq_batch(nib,[(start,end,strand)],mask)[0]


def get_nib_header(nib_fn,start=0,end=-1,strand='+',name=None,dbHeader=None,tbaHeader=None) :
    '''Method for constructing fasta headers compliant with nibFrag utility'''
    headers = get_nib_header_batch(nib,[(start,end,strand,name,dbHeader,tbaHeader),])
    return headers[0]


def get_nib_header_batch(nib,queries) :
    '''Batch method for creating nibFrag headers.  *queries* is a list of at most
    6-tuples (start,end,strand,name,dbHeader,tbaHeader) representing queries as
    specified by the original nibFrag utility.  Only start, end, and strand
    fields are required.'''

    nib_path, nib_f = _nib_fd(nib)

    nib_dir,nib_fn,nib_base,nib_ext = get_file_parts(nib_path)
    nbases = validate_nib_file(nib)
    headers = []
    header_tmpl = '>%(name)s%(db)s\n'

    for rec in queries :

        # set some defaults if they are not supplied
        rec = list(rec)
        rec.extend([None]*(6-len(rec)))
        start, end, strand, name, dbHeader, tbaHeader  = rec

        if end == -1 :
            end = nbases
        fields = {}
        fields['name'] = nib_path+':%d-%d'%(start,end) if not name else name
        fields['db'] = ''

        if tbaHeader :
            # ignored for some reason in nibFrag when tbaHeader supplied and dbHeader is not
            fields['name'] = '' if not dbHeader else fields['name']
            fields['db'] = '%s.%s:%d-%d of %d'%(tbaHeader,nib_base,start,end,nbases)
        if dbHeader :
            fields['db'] = ':%s.%s:%d-%d:%s:%d'%(dbHeader,nib_base,start,end,strand,nbases)

        headers.append(header_tmpl%fields)

    return headers


def validate_nib_file(nib) :
    '''Validate .nib file header, returning number of bases indicated if successful.
    *nib* argument is either a filename or an open file object.
    '''

    nib_fn, nib_f = _nib_fd(nib)

    # first 4 bytes are a nib file signature
    #TODO - consider attempting to figure out byte order to make truly cross platform
    def_sig = 0x6BE93D3A
    sig = struct.unpack('=l',nib_f.read(4))[0]
    if def_sig != sig :
        raise NibException('Invalid nib file signature in %s, found %s, expected \
        %s, perhaps .nib file as not created on this platform?\n\nnibFrag style \
        error: %s is not not a good .nib file.'%(nib_fn,hex(sig),hex(def_sig),nib_fn))

    # second 4 bytes are number of bases in sequence
    nbases = struct.unpack('=l',nib_f.read(4))[0]

    return nbases


def get_nib_seq_batch(nib,queries,mask=NOMASK) :
    '''Extract subsequence from .nib file like Jim Kent's nibFrag utility.

    Extract the nucleotide substrings defined by the closed intervals in *queries*
    from the sequence found in *nib*.  *nib* argument is either a filename or
    an open file object.  Entries in *queries* are 3-tuples defining (start,end,strand)
    sequence coordinates. Sequences are returned in order in a list as
    strings.  *mask* parameter has the following possible values:

    chipsequtil.nib.NOMASK -- masked positions are not indicated (default)
    chipsequtil.nib.MASK -- masked positions are capitalized, normal bases lower case
    chipsequtil.nib.NOMASK -- masked positions are replaced with Ns
    '''

    nib_fn, nib_f = _nib_fd(nib)

    nbases = validate_nib_file(nib_f)

    # rest of file is sequence, with each nibble (4 bytes) being a base as \
    # follows (from http://genome.ucsc.edu/FAQ/FAQformat.html#format8) :
    #
    # 0 - T
    # 1 - C
    # 2 - A
    # 3 - G
    # 4 - N
    #
    # The most significant bit in a nibble is set if the base is masked
    trans_nuc = 'tcagn'

    # start translating the nibbles into nucleotides
    def trans_nib(nib) :
        nuc = trans_nuc[nib&7]
        mask_bit = nib & 8
        if mask in [MASK,HARDMASK] and mask_bit == 0 :
            return nuc.upper()
        if mask == HARDMASK and mask_bit == 1 :
            return 'N'
        return nuc

    headers = [] # stores headers
    seqs = [] # stores sequences

    # sort the coords so we can walk most efficiently through the file
    queries.sort()

    for start, end, strand in queries :

        if start < 0 :
            raise NibException('Received negative start coordinate, this may '\
                               'indicate a region on mitochondrial DNA that '\
                               'spans reference sequence start and end.  This '\
                               'utility cannot handle these cases, aborting. '\
                               'Requested interval: %s (%d,%d)'%(nib_fn,start,end))

        start, end = map(int,(start,end))

        # end == -1 means caller wants entire sequence
        if end == -1  :
            end = nbases

        if any([nbases < c for c in [start,end]]) :
            raise NibException(('Requested slice (%(start)d,%(end)d) not compatible ' \
            'with sequence of length %(nbases)d in %(nib_fn)s, aborting\n\nnibFrag '\
            'style error: nib read past end of file (%(start)d %(end)d) in file: '\
            '%(nib_fn)s')%{'start':start,'end':end,'nbases':nbases,'nib_fn':nib_fn})

        # figure out how many bytes to read through
        start_byte,rem_byte = start/2,start%2

        # calculate where we need to move to in the file from the current location
        # + 8 is from the 2*4 bytes header info in the .nib format
        byte_offset = start_byte-nib_f.tell() + 8
        nib_f.seek(byte_offset,1) # seek forward to the beginning byte from current location
        seq_bytes,seq_rem_byte = int(math.ceil((end-start+rem_byte)/2.)),(end+1)%2
        seq_bytes = nib_f.read(seq_bytes+seq_rem_byte)

        # start translating the bytes
        seq = StringIO() # we use StringIO because it is more efficient than concatenating strings
        for c in seq_bytes :
            c_byte = struct.unpack('=b',c)[0]

            # higher nibble
            c_nib = (c_byte & (15<<4))>>4
            nuc = trans_nib(c_nib)
            seq.write(nuc)

            # lower nibble
            c_nib = int(c_byte) & 15
            nuc = trans_nib(c_nib)
            seq.write(nuc)

        # final nucleotide sequence
        seq_str = seq.getvalue()

        # if we're reading to the end, don't clip anything
        if end != nbases :
            # if the coordinate requested was not on a byte boundary, adjust
            if rem_byte == 1 :
                seq_str = seq_str[1:]
            if seq_rem_byte == 1 :
                seq_str = seq_str[:-1]

            # nibFrag apparently uses zero-based indexing, clip off one base
            seq_str = seq_str[:-1]
        seq.close()

        # adjust strand
        if strand == '-' :
            seq_str = reverse_complement(seq_str)
        seqs.append(seq_str)

    return seqs


class SeqDBException(Exception): pass
class NibDBException(Exception): pass


class SeqDB(object) :
    '''Base class for different kinds of sequence databases.  Does nothing,
    implement subclasses.  Constructor rovides _db_map and db_info class members.'''
    def __init__(self) :
        self._db_map = {}
        self.db_info = dd(dict)

    def get_seq(self,*args, **kwargs) :
        raise SeqDBException('Base class SeqDB has no get_seq implementation')


class NibDB(SeqDB) :
    '''Class providing an interface to a set of .nib files as created by faToNib
    in Jim Kent's software suite.

    Sequences are identified by the basename of the .nib file without the .nib
    extension, e.g. chr1.nib is identified as chr1.

    Some potentially useful information about the entries in the database is
    stored in the *nib_info* dictionary.
    '''

    def __init__(self,nib_fns=[],nib_dirs=[]) :
        '''*nib_fns* is a list of paths to specific .nib files desired for the
        NibDB.  *nib_dirs* is a list of paths to directories containing .nib
        files such that every .nib file in the directories is added to the NibDB.
        Explicitly passed files take precedence over those found in directories
        when sequence names collide.
        '''
        SeqDB.__init__(self)

        # find all *.nib files in the directories passed
        if isinstance(nib_dirs,str) : # user just provided single directory
            nib_dirs = [nib_dirs]

        dir_nibs = []
        for d in nib_dirs :
            dir_nibs.extend(glob.glob(os.path.join(d,'*.nib')))

        if isinstance(nib_fns,str) :
            nib_fns = [nib_fns]
        # for each .nib found, add to db
        # if there is a collision of names, those specified in files (not dirs)
        # takes precedence without warning
        for fn in dir_nibs+nib_fns :

            # open the nib file
            nib_path,nib_fn,nib_base,nib_ext = get_file_parts(fn)
            fn, nib_f = _nib_fd(fn)
            self._db_map[nib_base] = nib_f

            # store some info
            self.db_info[nib_base]['path'] = fn
            nbases = validate_nib_file(self._db_map[nib_base])
            self.db_info[nib_base]['nbases'] = nbases

    def __del__(self) :
        '''import this
        ...Explicit is better than implicit...
        '''
        for nib_f in self._db_map.values() :
            nib_f.close()

    def _get_db_map(self,name) :
        '''Gets appropriate file handle for the requested name, raises NibDBException
        if it cannot be found'''
        try :
            return self._db_map[name]
        except KeyError :
            raise NibDBException('Sequence name %s not found in NibDB'%name)

    def get_fasta(self,name,start=0,end=-1,strand='+',mask=NOMASK) :
        '''Get the fasta record for the specified arguments, returns (header,sequence)
        tuple.'''

        nib_f = self._get_db_map(name)
        return get_nib(nib_f,start,end,strand,mask)

    def get_fasta_batch(self,recs,mask=NOMASK) :
        '''Batch version of *get_fasta* method.  *recs* is a list of lists/tuples
        with (<chromo>,<start>,<end>,<strand>). Returns list of (header,sequence)
        tuples in the same sequence as the input records.'''

        # gather the records for each chromosome together
        chrom_recs = dd(list)
        for i,r in enumerate(recs) :
            chrom_recs[r[0]].append((i,r)) # recs are (index,<tuple>)

        # extract sequences
        all_chrom_recs = []
        for chrom, rec_list in chrom_recs.items() :
            # sorted lists make sequence extraction efficient
            rec_list.sort(key=lambda x: x[1][1]) # recs are (index,<tuple>)

            # separate indexes from records, extract for this chromo
            indexes, c_recs = zip(*rec_list)

            # get_nib_batch requires list of (<start>,<end>,<strand>) tuples, remove
            # chromo in first position
            c_recs = [r[1:] for r in c_recs]

            nib_f = self._get_db_map(chrom)
            headers, seqs = get_nib_batch(nib_f,c_recs,mask)

            # return the sequences to a (index,(header,sequence)) list
            all_chrom_recs.extend(zip(indexes,zip(headers,seqs)))

        # put the sequences back in the original order
        all_chrom_recs.sort(key=lambda x: x[0]) # recs are (index,<tuple>) again
        indexes, recs = zip(*all_chrom_recs)

        return zip(*recs)

    def get_fasta_from_bed(self,bed,mask=NOMASK) :
        '''Accepts either a chipsequtil.BEDFile instance or a filename for a BED
        file (used to construct a BEDFile instance) and returns the fasta
        records for all records in order.'''

        # determine if *bed* is a filename or a BEDFile
        if isinstance(bed,str) : # filename
            bed = BEDFile(bed)

        # construct the records
        recs = []
        for rec in bed :
            if rec['chrom'].lower().startswith('track') : # track line, skip
                continue
            recs.append((rec['chrom'],int(rec['chromStart']),int(rec['chromEnd']),rec['strand']))

        return self.get_fasta_batch(recs,mask)

    def get_seq(self,name,start=0,end=-1,strand='+',mask=NOMASK) :
        '''Extract sequence from sequence *name*. Other arguments are passed
        directly to *get_nib_seq* function.'''

        nib_f = self._get_db_map(name)
        return get_nib_seq(nib_f,start,end,strand,mask)
