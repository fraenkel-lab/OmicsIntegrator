'''
Fasta.py -- Very efficient code for loading biological sequences in Fasta format
            python dictionaries.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
'''

import sys, re, os, random
from   gzip import GzipFile

def keys(filename,key_func=None):
    '''
    Fasta.keys(filename,key_func=None)
    ----------------------------------
    Return the ids in a Fasta file.  Same as Fasta.ids(file,key_func=None)
    key_func is a function (or lambda expression that extracts the key.
    Default is to take first word separated by whitespace.

    For example:

    key_func=lambda x: x.split('|')[3]

    Would use the 4th token separated by the "|" symbol
    '''
    FID = open(filename)
    first = FID.readline()[0]
    FID.seek(0)
    if not re.search('\.fa|\.fsa|\.fasta',filename) or (first != '>'): #Fasta file?
        ids = [x.strip().split()[0] for x in FID.readlines() if (x.strip() and x[0] != '#') ]
        FID.close()
        return ids
    return file2dict(filename,key_func=key_func).keys()

def ids(filename,key_func=None):
    '''
    Fasta.ids(filename,key_func=None)
    ---------------------------------
    Return the ids in a Fasta file.  Same as Fasta.keys(file,key_func=None)
    key_func is a function (or lambda expression that extracts the key.
    Default is to take first word separated by whitespace.

    For example:

    key_func=lambda x: x.split('|')[3]

    Would use the 4th token separated by the "|" symbol
    '''
    return keys(filename,key_func=key_func)

def seqs(filename):
    '''
    Fasta.seqs(filename)
    --------------------
    Return a list of the sequences contained in the file.
    '''
    return load(filename).values()

def load(filename,key_func=None):
    '''
    Fasta.load(filename,key_func=None)
    ---------------------------------
    Load the file "filename" as a dictionary of sequences, indedex according
    to key_func. Default is to take first word separated by whitespace.

    For example:

    key_func=lambda x: x.split('|')[3]

    Would use the 4th token separated by the "|" symbol
    '''
    return fasta2dict(filename,key_func=key_func)

def file2dict(filename,key_func=None):
    '''
    Fasta.file2dict(filename,key_func=None)
    --------------------------------------
    Synonymous with Fasta.load().  See documentation for Fasta.load().
    '''
    return fasta2dict(filename,'WANT_DICT',key_func)

def fasta2dict(filename, want_dict = 'YES',key_func=None):
    '''
    Fasta.fasta2dict(filename, want_dict = 'YES',key_func=None)
    ----------------------------------------------------------
    Very fast Fasta Loader.  Used internally.  You should be using
    Fasta.load() or Fasta.seqs() instead.
    '''
    D = {}
    if filename[-3:] == '.gz': FH = GzipFile(filename)
    else:                      FH = open(filename,'r')
    chunks = FH.read().split('>')
    for chunk in chunks[1:]:
        lines  = chunk.split('\n')
        raw_id = lines[0]
        seq    = ''.join(lines[1:])
        try:
            if not key_func:  key = raw_id.split()[0]
            else:            key = key_func(raw_id)
        except:
            print raw_id
            sys.stdout.flush()
            sys.exit(1)
        D[key] = seq
    if want_dict: return D
    else:         return D.values()


def delN(fsaD):
    '''
    Fasta.delN(fsaD)
    ----------------
    Remove any entries in the Fasta-derived dictionary that have any DNA
    ambiguity codes within.  Reports ids of deleted sequences.
    '''
    for key,seq in fsaD.items():
        seq = re.sub("^N*","",seq)
        seq = re.sub("N*$","",seq)
        if re.search('[NRYKMSWBDHV]',seq):
            print 'deleting ',key
            del fsaD[key]
        else:
            fsaD[key] = seq

def find(name,pathhint=None):
    '''
    Fasta.find(name,pathhint=None)
    ------------------------------
    Find a ".fsa" file with a similar name to the supplied file.

    For example, given "GAL4_YPD.meme," this function will look in
    the current directory, then the parent directory or the
    optinal "hint" directory for a file with the name "GAL4_YPD.fsa"
    
    '''
    exists = os.path.exists
    root   = re.sub('\.\w*$','',name)
    smroot = re.sub('_.$'   ,'',root)
    if pathhint: pathhint = pathhint + '/'
    tail   = root.split('/')[-1]
    parent = '/'.join(name.split('/')[:-2])
    if parent: parent = parent + '/'

    if re.search('\.fsa$',name):
        if exists(name):
            return name
        elif pathhint and exists(pathhint + tail + '.fsa'):
            return pathhint + tail + '.fsa'
    else:
        if exists(root + '.fsa'):
            return root + '.fsa'
        elif pathhint and exists(pathhint + tail + '.fsa'):
            return pathhint + tail + '.fsa'
        elif exists(parent + tail + '.fsa'):
            return parent + tail + '.fsa'
        elif exists(smroot + '.fsa'):
            return smroot + '.fsa'
        elif pathhint and exists(pathhint + smroot + '.fsa'):
            return pathhint + smroot + '.fsa'
    print '## ! Could not find fsa file for %s'%name
    return None

def write(D,filename,linelen=70):
    '''
    Fasta.write(D,filename,linelen=70)
    ----------------------------------
    Write dictionary of sequences out to a file.  Optional
    linelen argument specifies how many sequence characters
    are allowed on each line.
    '''
    F = open(filename,'w')
    F.write(text(D,linelen=linelen))
    F.close()


def text(D,toupper=0,linelen=70):
    '''
    Fasta.text(D,toupper=0,linelen=70)
    ----------------------------------
    Utility function for generating Fasta-formatted output from a dictonary
    of sequences.  toupper specifies if all sequences should be capitalized,
    and linelen specifies how many sequence characters are allowed on each line.

    Returns a single string.
    '''
    sA = []
    keys = D.keys()
    keys.sort()
    for id in keys:
        if toupper: seq = D[id].upper()
        else:       seq = D[id]
        sA.append(">%s\n"%id)
        for i in range(0,len(seq),linelen):
            sA.append(seq[i:i+linelen] + '\n')
    s = ''.join(sA)
    return(s[:-1])


def random_subset(filename_or_seqD,target_count=30):
    '''
    Fasta.random_subset(filename_or_seqD,target_count=30)
    -----------------------------------------------------
    Pick a subset of entries at random from the specificied
    input, and return a dictionary.  The input may be either
    a filename or a dictionary of sequences.

    "target_count" is the desired size of the subset.
    '''
    
    if type(filename_or_seqD) == type({}):
        seqD = filename_or_seqD
    else:
        seqD = file2dict(filename_or_seqD)
    ids      = seqD.keys()
    newD     = {}
    count    = 0
    numseqs  = len(ids)
    while count < target_count:
        randomid = ids[int(random.random()*numseqs)]
        if newD.has_key(randomid): continue
        newD[randomid] = seqD[randomid]
        count = count + 1
    return newD
    

def random_split(filename_or_seqD,frac=0.5):
    '''
    Fasta.random_split(filename_or_seqD,frac=0.5)
    --------------------------------------------
    Randomly partition a fasta-derived dictioary. The input may be either
    a filename or a dictionary of sequences.  The "frac" argument
    specifies the ratio of number of sequences.

    Returns two dictionaries.
    '''
    if type(filename_or_seqD) == type({}):
        seqD = filename_or_seqD
    else:
        seqD    = file2dict(filename_or_seqD)
    newD    = {}
    remainD = {}
    ids = seqD.keys()
    targetcount = int(frac * len(ids))
    count = 0

    while count< targetcount:
        randomid = ids[int(random.random()*len(ids))]
        if newD.has_key(randomid): continue
        newD[randomid] = seqD[randomid]
        count = count + 1

    for id in seqD.keys():
        if newD.has_key(id): continue
        remainD[id] = seqD[id]

    return newD, remainD
