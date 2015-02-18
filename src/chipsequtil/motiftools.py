"""
There is a large number of functions and member fucntions here.  To get started,
a motif can be instantiated by providing an ambiguity code, a set of aligned DNA
sequences, or from matrices of counts, probabilities or log-likelihoods (akaPSSMs).

>>> m =  MotifTools.Motif_from_text('TGAAACANNSYWT')
>>> print m.oneletter()
TGAAACA..sywT 

Lower case reflects lower information content.  For a more detailed view of the distribution
of information, try this::

    >>> m.textlogo()
    #                -- 2.30 bits
    #
    #  TGAAACA     T
    #  TGAAACA     T
    #  TGAAACA     T
    #  TGAAACA     T
    #  TGAAACA  CCAT
    #  TGAAACA  CCAT
    #  TGAAACA  GTTT
    #  TGAAACA  GTTT -- 0.23 bits
    #  -------------
    #  TGAAACA..sywT


Motif objects may be manipulated largely like text strings (with pythonic
indexing)::

    >>> print m[4:5].oneletter
    A 
    >>> print m[4:7].oneletter
    ACA 
    >>> print (m[4:7] + m[1:2]).oneletter
    ACAG
    >>> print (m[4:7] + m[1:7]).oneletter
    ACAGAAACA

and even padded with blanks::

    >>> print  m[-4:7]
    ...TGAAACA

.. Copyright (2005) Whitehead Institute for Biomedical Research
.. All Rights Reserved

Author: David Benjamin Gordon

Modified by: Adam Labadorf

"""
import copy
import math
import os
import pickle
import re
import string
import sys
import tempfile

pysum = sum

from random import random,shuffle
from subprocess import call

from chipsequtil import reverse_complement
class MotifToolsException(Exception) : pass

one2two = {  'W':'AT',    'M':'AC',   'R':'AG',
             'S':'CG',    'Y':'CT',   'K':'GT'}
two2one = { 'AT': 'W',   'AC': 'M',  'AG': 'R',
            'CG': 'S',   'CT': 'Y',  'GT': 'K'}
revcomp = { 'A':'T',      'T':'A',    'C':'G',   'G':'C',
            'W':'W',      'S':'S',    'K':'M',   'M':'K',
            'Y':'R',      'R':'Y',    'N':'N',
            'B':'N', 'D':'N', 'H':'N', 'V':'N', ' ':'N'}  #[12-11-02] Needs fixing

ACGT = list('ACGT')
YEAST_BG = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast default background freqs

revcomplement_memo = {'A':'T'}
revcompTBL = string.maketrans("AGCTagctWSKMYRnN", "TCGAtcgaWSMKTYnN")
def revcomplement(seq):
    """A quick reverse-complement routine that memo-izes queries, understands
    IUPAC ambiguity codes, and preserves case."""
    global revcomplement_memo
    try:
        rc = revcomplement_memo[seq]
    except KeyError:
        #_t = map(lambda x,D=revcomp: D[x], seq)
        #get = revcomp.get
        #_t = map(get, seq)
        _t = list(seq.translate(revcompTBL))
        _t.reverse()
        rc = ''.join(_t)
        revcomplement_memo[seq] = rc
        revcomplement_memo[rc]  = seq
    return rc


def Motif_from_ll(ll):
    """Constructs a motif object from a log-likelihood matrix, which is in the
    form of a list of dictionaries."""
    m = Motif(None,None)
    m.compute_from_ll(ll)
    return m

def Motif_from_counts(countmat,beta=0.01,bg={'A':.25,'C':.25,'G':.25,'T':.25}):
    """
    Construct a Motif object from a matrix of counts (or probabilities or frequencies).
    A default set of uniform background frequencies may be overridden.

    beta refers to the number of pseudocounts that should be distributed over each position
    of the PSSM."""
    m = Motif('',bg)
    m.compute_from_counts(countmat,beta)
    return m

def Motif_from_text(text,beta=0.05,source='',bg=None):
    """Construct a Motif object from a text string constructed from IUPAC
    ambiguity codes. 

    A default set of uniform background frequencies may be overridden with
    a dictionary of the form {'A':.25,'C':.25,'G':.25,'T':.25}).

    beta refers to the number of pseudocounts that should be distributed over each position
    of the PSSM."""
    if not bg: bg={'A':.25,'C':.25,'G':.25,'T':.25}
    m = Motif('',bg)
    m.compute_from_text(text,beta)
    m.source = source
    return m

def copy(motif):
    """Utility routine for copying motifs"""
    a = copy.deepcopy(motif)
    #a.__dict__ = motif.__dict__.copy()
    return a

class Motif:
    """A pssm model, with scanning, storing, loading, and other operations. A
    uniform nucleotide background is assumed if none is provided."""
    def __init__(self,list_of_seqs_or_text=[],backgroundD=None):
        self.MAP       = 0
        self.evalue    = None
        self.oneletter = ''
        self.nseqs     = 0
        self.counts    = []
        self.width     = 0
        self.fracs     = []
        self.logP      = []
        self.ll        = []
        self.bits      = []
        self.totalbits = 0
        self.maxscore  = 0
        self.minscore  = 0
        self.pvalue      = 1
        self.pvalue_rank = 1
        self.church      = None
        self.church_rank = 1
        self.Cpvalue     = 1
        self.Cpvalue_rank= 1
        self.Cchurch     = 1
        self.Cchurch_rank= 1
        self.binomial    = None
        self.binomial_rank=1
        self.E_seq       = None
        self.frac        = None
        self.E_site      = None
        self.E_chi2      = None
        self.kellis      = None
        self.MNCP        = None
        self.ROC_auc     = None
        self.realpvalue  = None
        self.Cfrac       = None
        self.CRA         = None
        self.valid     = None
        self.seeddist  = 0
        self.seednum   = -1
        self.seedtxt   = None
        self.family    = None
        self.source    = None
        self.threshold = None
        self._bestseqs = None
        self.bgscale   = 1
        self.best_pvalue = None
        self.best_factor = None
        self.gamma     = None
        self.nbound    = 0
        self.matchids  = []
        self.overlap   = None
        self.cumP      = []
        self.numbound      = 0
        self.nummotif      = 0
        self.numboundmotif = 0
        self.dataset = None
        self.bgfile = None
        self.cverror = None
        self.beta = None
        self.match_thresh = None
        self.progscore = None
        if backgroundD:
            self.background = backgroundD
        else:
            #self.background = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast Default
            self.background = {'A':.25,'C':.25,'G':.25,'T':.25} # uniform background

        if type(list_of_seqs_or_text) == type(''):
            self.seqs = []
            text = list_of_seqs_or_text
            self.compute_from_text(text)
        else:
            self.seqs = list_of_seqs_or_text
        if self.seqs:
            self._parse_seqs(list_of_seqs_or_text)
            self._compute_ll()
            self._compute_oneletter()
            #self._compute_threshold(2.0)

    def __repr__(self):
        return "%s (%d)"%(self.oneletter, self.nseqs)

    def __str__(self):
        return "%s (%d)"%(self.oneletter, self.nseqs)

    def summary(self):
        """return a text string one-line summary of motif and its metrics"""
        m = self
        txt = "%-34s (Bits: %5.2f  MAP: %7.2f   D: %5.3f  %3d)  E: %7.3f"%(
            m, m.totalbits, m.MAP, m.seeddist, m.seednum, nlog10(m.pvalue))
        if m.binomial!=None:  txt = txt + '  Bi: %6.2f'%(nlog10(m.binomial))
        if m.church != None:  txt = txt + '  ch: %6.2f'%(nlog10(m.church))
        if m.frac   != None:  txt = txt + '  f: %5.3f'%(m.frac)
        if m.E_site != None:  txt = txt + '  Es: %6.2f'%(nlog10(m.E_site))
        if m.E_seq  != None:  txt = txt + '  Eq: %6.2f'%(nlog10(m.E_seq))
        if m.MNCP   != None:  txt = txt + '  mn: %6.2f'%(m.MNCP)
        if m.ROC_auc!= None:  txt = txt + '  Ra: %6.4f'%(m.ROC_auc)
        if m.E_chi2 != None:
            if m.E_chi2 == 0: m.E_chi2=1e-20
            txt = txt + ' x2: %5.2f'%(nlog10(m.E_chi2))
        if m.CRA    != None:  txt = txt + '  cR: %6.4f'%(m.CRA)
        if m.Cfrac  != None:  txt = txt + '  Cf: %5.3f'%(m.Cfrac)
        if m.realpvalue != None: txt = txt + '  P: %6.4e'%(m.realpvalue)
        if m.kellis != None:  txt = txt +  '  k: %6.2f'%(m.kellis)
        if m.numbound      :  txt = txt +  '  b: %3d'%(m.numbound)
        if m.nummotif      :  txt = txt +  '  nG: %3d'%(m.nummotif)
        if m.numboundmotif :  txt = txt +  '  bn: %3d'%(m.numboundmotif)

        return txt

    def minimal_raw_seqs(self):
        '''return minimal list of seqs that represent consensus '''
        seqs = [[], []]
        for letter in self.oneletter:
            if one2two.has_key(letter):
                seqs[0].append(one2two[letter][0])
                seqs[1].append(one2two[letter][1])
            else:
                seqs[0].append(letter)
                seqs[1].append(letter)
        if ''.join(seqs[0]) == ''.join(seqs[1]):
            return  [''.join(seqs[0])] 
        else:
            return  [''.join(seqs[0]), ''.join(seqs[0])] 
    def _compute_oneletter(self):
        """set the oneletter member variable"""
        letters = []
        for i in range(self.width):
            downcase = None
            if self.bits[i] < 0.25:
                letters.append('.')
                continue
            if self.bits[i] < 1.0: downcase = 'True'
            tups = [(self.ll[i][x],x) for x in ACGT if self.ll[i][x] > 0.0]
            if not tups:  #Kludge if all values are negative (can this really happen?)
                tups = [(self.ll[i][x],x) for x in ACGT]
                tups.sort()
                tups.reverse()
                tups = [tups[0]]
                downcase = 'True'
            tups.sort()      #Rank by LL
            tups.reverse()
            bases = [x[1] for x in tups[0:2]]
            bases.sort()
            if len(bases) == 2: L = two2one[''.join(bases)]
            else:               L = bases[0]
            if downcase: L = L.lower()
            letters.append(L)
        self.oneletter = ''.join(letters)
    def _parse_seqs(self, LOS):
        """build a matrix of counts from a list of sequences"""
        self.nseqs = len(LOS)
        self.width = len(LOS[0])
        for i in range(self.width):
            Dc = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}
            for seq in LOS:
                key = seq[i]
                Dc[key] = Dc[key] + 1
            del(Dc['N'])
            self.counts.append(Dc)

    def _compute_ll(self):
        """compute the log-likelihood matrix from the count matrix"""
        self.fracs = []
        self.logP  = []
        self.ll    = []
        for i in range(self.width):

            Dll  = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            Df   = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            DlogP= {'A': 0, 'C': 0, 'T': 0, 'G': 0}

            for nuc in self.counts[i].keys():

                #print i,nuc,self.counts[i][nuc],self.nseqs
                # Dll[nuc] = log2( position nucleotide count/background sequence count )
                # Dll[nuc] = log2( (count[nuc]+bgscale*bg[nuc])/(bg[nuc]*(num_seqs+bgscale)) )

                pos_nuc_count = self.counts[i][nuc] + self.bgscale*self.background.get(nuc,0.)
                adj_all_nuc_count = (self.nseqs + self.bgscale) * self.background.get(nuc,1e-10)

                Dll[nuc] = math.log(pos_nuc_count/adj_all_nuc_count,2)

                Pij = self.counts[i][nuc] / float(self.nseqs)
                Df [nuc] = Pij
                if Pij > 0:
                    DlogP[nuc]  = math.log(Pij) / math.log(2.)
                else:
                    DlogP[nuc]  = -100  #Near zero

            self.fracs.append(Df)
            self.logP.append (DlogP)
            self.ll.append   (Dll)
        self.P = self.fracs
        self._compute_bits()
        self._compute_ambig_ll()
        self._maxscore()


    def compute_from_ll(self,ll):
        """build motif from an inputed log-likelihood matrix

        (This function reverse-calculates the probability matrix and background frequencies
        that were used to construct the log-likelihood matrix)
        """
        self.ll    = ll
        self.width = len(ll)
        self._compute_bg_from_ll()
        self._compute_logP_from_ll()
        self._compute_ambig_ll()
        self._compute_bits()
        self._compute_oneletter()
        self._maxscore()

    def _computeP(self):
        """compute the probability matrix (from the internal log-probability matrix)"""
        P = []
        for i in range(self.width):
            #print i,
            _p = {}
            for L in ACGT: _p[L] = math.pow(2.0,self.logP[i][L])
            P.append(_p)
        #print
        self.P = P

    def _compute_bits(self):
        """set m.totbits to the number of bits and m.bits to a list of bits at
        each position"""
        bits = []
        totbits = 0
        bgbits  = 0
        bg      = self.background
        UNCERT  = lambda x: x*math.log(x)/math.log(2.0)
        for letter in ACGT:
            bgbits = bgbits + UNCERT(bg[letter])
        for i in range(self.width):
            tot = 0
            for letter in ACGT:
                Pij = pow(2.0, self.logP[i][letter])
                tot = tot + UNCERT(Pij)
                #bit = Pij * self.ll[i][letter]
                #if bit > 0:
                #    tot = tot + bit
            #print tot, bgbits, tot-bgbits
            bits.append(max(0,tot-bgbits))
            totbits = totbits + max(0,tot-bgbits)
        self.bits = bits
        self.totalbits = totbits

        
    def denoise(self,bitthresh=0.5):
        """set low-information positions (below bitthresh) to Ns"""
        for i in range(self.width):
            tot = 0
            for letter in ACGT:
                if self.logP:
                    Pij = pow(2.0, self.logP[i][letter])
                else:
                    Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                if Pij > 0.01:
                    bit = Pij * self.ll[i][letter]
                    tot = tot + bit
            if tot < bitthresh:  #Zero Column
                for letter in ACGT:
                    self.ll[i][letter] = 0.0
        self.compute_from_ll(self.ll)

    def giflogo(self,id,title=None,scale=0.8,info_str=''):
        """make a gif sequence logo"""
        return giflogo(self,id,title,scale)

    def printlogo(self,norm=2.3, height=10.0):
        """print a text-rendering of the Motif Logo

        norm
            maximum number of bits to show
        height
            number of lines of text to use to render logo
        """
        self._print_bits(norm,height)
    def print_textlogo(self,norm=2.3, height=8.0):
        """print a text-rendering of the Motif Logo

        norm
            maximum number of bits to show
        height
            number of lines of text to use to render logo
        """
        self._print_bits(norm,height)
    def _print_bits(self,norm=2.3, height=8.0):
        """print a text-rendering of the Motif Logo

        norm
            maximum number of bits to show
        height
            number of lines of text to use to render logo
        """
        bits   = []
        tots   = []
        str    = []
        for i in range(self.width):
            D = {}
            tot = 0
            for letter in ['A', 'C', 'T', 'G']:
                if self.logP:
                    Pij = pow(2.0, self.logP[i][letter])
                else:
                    Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                if Pij > 0.01:
                    '''Old'''
                    D[letter] = Pij * self.ll[i][letter]
                    #'''new'''
                    #Q = self.background[letter]
                    #D[letter] = ( Pij * math.log(Pij) - Pij * math.log(Q) ) / math.log(2.0)
                    '''for both old and new'''
                    tot = tot + D[letter]
            bits.append(D)
            tots.append(tot)
        for i in range(self.width):
            s = []
            _l = bits[i].keys()
            _l.sort(lambda x,y,D=bits[i]: cmp(D[y],D[x]))
            for key in _l:
                for j in range(int(bits[i][key] / norm * height)):
                    s.append(key)
            str.append(''.join(s))
        fmt = '%%%ds'%height
        print '#  %s'%('-'*self.width)
        for h in range(int(height)):
            sys.stdout.write("#  ")
            for i in range(self.width):
                sys.stdout.write((fmt%str[i])[h])
            if h == 0:
                sys.stdout.write(' -- %4.2f bits\n'%norm)
            elif h == height-1:
                sys.stdout.write(' -- %4.2f bits\n'%(norm/height))
            else:
                sys.stdout.write('\n')
        print '#  %s'%('-'*self.width)
        print '#  %s'%self.oneletter

    def _compute_ambig_ll(self):
        """extend log-likelihood matrix to include ambiguity codes
        e.g.  What the score of a 'S'?  Here we use the max of C and G."""
        for Dll in self.ll:
            for L in one2two.keys():
                Dll[L] = max(Dll[one2two[L][0]],  Dll[one2two[L][1]] )
            Dll['N'] = 0.0
            Dll['B'] = 0.0

    def compute_from_nmer(self,nmer,beta=0.001):  #For reverse compatibility
        """See compute_from_text.  Here for reverse compatibility"""
        self.compute_from_text(nmer,beta)

    def compute_from_text(self,text,beta=0.001):
        """compute a matrix values from a text string of ambiguity codes.
        Use Motif_from_text utility instead to build motifs on the fly."""
        prevlett = {'B':'A', 'D':'C', 'V':'T', 'H':'G'}
        countmat = []
        text = re.sub('[\.\-]','N',text.upper())
        for i in range(len(text)):
            D = {'A': 0, 'C': 0, 'T':0, 'G':0}
            letter = text[i]
            if letter in ['B', 'D', 'V', 'H']:  #B == no "A", etc...
                _omit = prevlett[letter]
                for L in ACGT:
                    if L != _omit: D[L] = 0.3333
            elif one2two.has_key(letter):  #Covers WSMYRK
                for L in list(one2two[letter]):
                    D[L] = 0.5
            elif letter == 'N':
                for L in D.keys():
                    D[L] = self.background[L]
            elif letter == '@':
                for L in D.keys():
                    D[L] = self.background[L]-(0.0001)
                D['A'] = D['A'] + 0.0004
            else:
                D[letter] = 1.0
            countmat.append(D)
        self.compute_from_counts(countmat,beta)

    def new_bg(self,bg):
        """change the ACGT background frequencies to those in the supplied dictionary.
        Recompute log-likelihood, etc. with new background.
        """
        counts = []
        for pos in self.logP:
            D = {}
            for L,lp in pos.items():
                D[L] = math.pow(2.0,lp)
            counts.append(D)
        self.background = bg
        self.compute_from_counts(counts,0)

    def addpseudocounts(self,beta=0):
        """add pseudocounts uniformly across the matrix"""
        self.compute_from_counts(self.counts,beta)

    def compute_from_counts(self,countmat,beta=0):
        """build a motif object from a matrix of letter counts."""
        self.counts  = countmat
        self.width   = len(countmat)
        self.bgscale = 0

        maxcount = 0
        #Determine Biggest column
        for col in countmat:
            tot = pysum(col.values())
            if tot > maxcount :
                maxcount = tot

        #Pad counts of remaining columns
        for col in countmat:
            tot = pysum(col.values())
            pad = maxcount - tot
            for L in col.keys():
                col[L] = col[L] + pad * self.background.get(L,0.)

        self.nseqs = maxcount
        nseqs = maxcount

        #Add pseudocounts
        if beta > 0:  
            multfactor = {}
            bgprob = self.background
            pcounts= {}
            for L in bgprob.keys():
                pcounts[L] = beta*bgprob[L]*nseqs 
            for i in range(self.width):
                for L in countmat[i].keys():
                    _t = (countmat[i][L] + pcounts[L]) #Add pseudo
                    _t = _t / (1.0 + beta)    #Renomalize
                    countmat[i][L] = _t

        #Build Motif
        self.counts = countmat
        self._compute_ll()
        self._compute_oneletter()
        self._maxscore()


    def _compute_bg_from_ll(self):
        """compute background model from log-likelihood matrix
        by noting that:   pA  + pT  + pC  + pG  = 1
                  and     bgA + bgT + bgC + bgG = 1
                  and     bgA = bgT,   bgC = bgG
                  and so  bgA = 0.5 - bgC
                  and     pA  = lA * bgA,  etc for T, C, G
                  so...
                         (lA + lT)bgA + (lC + lG)bgC          =  1
                         (lA + lT)bgA + (lC + lG)(0.5 - bgA)  =  1
                         (lA + lT - lC - lG)bgA +(lC +lG)*0.5 =  1
                          bgA                                 =  {1 - 0.5(lC + lG)} / (lA + lT - lC - lG)
        + Gain accuracy by taking average of bgA over all positions of PSSM
        """
        
        pow = math.pow
        bgATtot = 0
        nocount = 0
        near0   = lambda x:(-0.01 < x and x < 0.01)
        for i in range(self.width):
            _D = self.ll[i]
            ATtot = pow(2,_D['A']) + pow(2,_D['T'])
            GCtot = pow(2,_D['C']) + pow(2,_D['G'])
            if near0(_D['A']) and near0(_D['T']) and near0(_D['G']) and near0(_D['C']):
                nocount = nocount + 1
                continue
            if near0(ATtot-GCtot):     #Kludge to deal with indeterminate case
                nocount = nocount + 1
                continue
            bgAT   = (1.0 - 0.5*GCtot)/(ATtot - GCtot)
            if (bgAT < 0.1) or (bgAT > 1.1):
                nocount = nocount + 1
                continue
            bgATtot = bgATtot + bgAT
        if nocount == self.width:  #Kludge to deal with different indeterminate case
            self.background = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
            return
        bgAT = bgATtot / (self.width - nocount)
        bgGC = 0.5 - bgAT
        self.background = {'A':bgAT, 'C':bgGC, 'G':bgGC, 'T':bgAT}            
        
    def _compute_logP_from_ll(self):
        """compute self's logP matrix from the self.ll (log-likelihood)"""
        log = math.log
        logP = []
        for i in range(self.width):
            D = {}
            for L in ACGT:
                ''' if   ll = log(p/b) then
                       2^ll = p/b
                  and    ll = log(p) - log(b)
                  so log(p) = ll + log(b)'''
                #Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                D[L] = self.ll[i][L] + log(self.background[L])/log(2.)
            logP.append(D)
        self.logP = logP

    def _print_ll(self):
        """print log-likelihood (scoring) matrix"""
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%self.ll[i][L],
            print
    def _print_p(self):
        """print probability (frequency) matrix"""
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%math.pow(2,self.logP[i][L]),
            print
    def _print_counts(self):
        """print count matrix"""
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%self.counts[i][L],
            print
        
    def _maxscore(self):
        """sets self.maxscore and self.minscore"""
        total = 0
        lowtot= 0
        for lli in self.ll:
            total = total + max(lli.values())
            lowtot= lowtot+ min(lli.values())
        self.maxscore = total
        self.minscore = lowtot

    def _compute_threshold(self,z=2.0):
        """for Motif objects assembled from a set of sequence,
        compute a self.threshold with a z-score based on the distribution
        of scores in among the original input sequences.
        """
        scoretally = []
        for seq in self.seqs:
            matches,endpoints,scores = self.scan(seq,-100)
            scoretally.append(scores[0])
        ave,std = avestd(scoretally)
        self.threshold = ave - z *std
        #print '#%s: threshold %5.2f = %5.2f - %4.1f * %5.2f'%(
        #    self, self.threshold, ave, z, std)

    def bestscanseq(self,seq):
        """return score,sequence of the best match to the motif in the supplied sequence"""
        matches,endpoints,scores = self.scan(seq,-100)
        t = zip(scores,matches)
        t.sort()
        bestseq   = t[-1][1]
        bestscore = t[-1][0]
        return bestscore, bestseq
    
    def bestscore(self,seq):
        """return the score of the best match to the motif in the supplied sequence"""
        return m.bestscan(seq)

    def bestscan(self,seq):
        """return the score of the best match to the motif in the supplied sequence"""
        matches,endpoints,scores = self.scan(seq,-100)
        if not scores: return -100
        scores.sort()
        best = scores[-1]
        return best

    def matchstartorient(self,seq, factor=0.7):
        """returns list of (start,orientation) coordinate pairs of matches to
        the motif in the supplied sequence.  Factor is multiplied by m.maxscore
        to get a match threshold.
        """
        ans = []
        txts,endpoints,scores = self.scan(seq,factor=factor)
        for txt, startstop in zip(txts,endpoints):
            start, stop = startstop
            rctxt  = reverse_complement(txt)
            orient = (self.bestscore(txt,1) >= self.bestscore(rctxt,1))
            ans.append((start,orient))
        return ans

    def scan(self, seq, threshold = '', factor=0.7):
        """
        Scan the sequence.  Returns three lists: matching sequences, endpoints,
        and scores.  The value of 'factor' is multiplied by m.maxscore to get a
        match threshold if none is supplied
        """
        if len(seq) < self.width:
            return self._scan_smaller(seq,threshold)
        else:
            return self._scan(seq,threshold,factor=factor)

    def scansum(self,seq,threshold = -1000):
        """
        Sum of scores over every window in the sequence.  Returns
        total, number of matches above threshold, average score, sum of exp(score)
        """
        ll = self.ll
        sum = 0
        width        = self.width
        width_r      = range(width)
        width_rcr    = range(width-1,-1,-1)
        width_ranges = zip(width_r,width_rcr)
        seqcomp      = seq.translate(revcompTBL)

        total = 0
        hits  = 0
        etotal= 0
        for offset in range(len(seq)-width+1):
            total_f = 0
            total_r = 0
            for i,ir in width_ranges:
                pos = offset+i
                total_f = total_f + ll[i][    seq[pos]]
                total_r = total_r + ll[i][seqcomp[pos]]
            total_max = max(total_f,total_r)
            if total_max >= threshold:
                total = total + total_max
                etotal = etotal + math.exp(total_max)
                hits  = hits + 1
            if not hits:
                ave = 0
            else:
                ave = float(total)/float(hits)
        return total,hits,ave,math.log(etotal)

    def score(self, seq, fwd='Y'):
        """returns the score of the first w-bases of the sequence, where w is the motif width."""
        matches, endpoints, scores = self._scan(seq,threshold=-100000,forw_only=fwd)
        return scores[0]

    def bestscore(self,seq, fwd=''):
        """returns the score of the best matching subsequence in seq."""
        matches, endpoints, scores = self._scan(seq,threshold=-100000,forw_only=fwd)
        if scores: return max(scores)
        else:      return -1000

    def _scan(self, seq,threshold='',forw_only='',factor=0.7):
        """internal tility function for performing sequence scans"""
        ll = self.ll #Shortcut for Log-likelihood matrix
        if not threshold: threshold = factor * self.maxscore
        
        #print '%5.3f'%(threshold/self.maxscore)
        matches       = []
        endpoints     = []
        scores        = []
        width         = self.width
        width_r       = range(width)
        width_rcr     = range(width-1,-1,-1)
        width_ranges  = zip(width_r,width_rcr)

        seqcomp = seq.translate(revcompTBL)

        for offset in range(len(seq)-self.width+1):    #Check if +/-1 needed
            total_f = 0
            total_r = 0
            for i,ir in width_ranges:
                pos = offset+i
                total_f = total_f + ll[i ][    seq[pos]]
                total_r = total_r + ll[ir][seqcomp[pos]]

            if 0 and total_f > 1:
                for i,ir in width_ranges:
                    print seq[offset+i],'%6.3f'%ll[i ][        seq[offset+i] ],'   ',
                print '= %7.3f'%total_f
                
            if 0:
                print "\t\t%s vs %s: F=%6.2f R=%6.2f %6.2f %4.2f"%(seq[offset:offset+self.width],
                                                                   self.oneletter,total_f,total_r,
                                                                   self.maxscore,
                                                                   max([total_f,total_r])/self.maxscore)
            if total_f > threshold and ((total_f > total_r) or forw_only):
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_f)
                matches.append(seq[offset:offset+self.width])
            elif total_r > threshold:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_r)
                matches.append(seq[offset:offset+self.width])
        return matches,endpoints,scores
    def _scan_smaller(self, seq, threshold=''):
        """internal utility function for performing sequence scans. The sequence
        is smaller than the PSSM.  Are there good matches to regions of the PSSM?"""
        ll = self.ll #Shortcut for Log-likelihood matrix
        matches   = []
        endpoints = []
        scores    = []
        w         = self.width
        for offset in range(self.width-len(seq)+1):    #Check if +/-1 needed
            maximum = 0
            for i in range(len(seq)):
                maximum = maximum + max(ll[i+offset].values())
            if not threshold: threshold = 0.8 * maximum
            total_f = 0
            total_r = 0
            for i in range(len(seq)):
                total_f = total_f + ll[i+offset      ][        seq[i] ]
                total_r = total_r + ll[w-(i+offset)-1][revcomp[seq[i]]]
            if 0:
                print "\t\t%s vs %s: F=%6.2f R=%6.2f %6.2f %4.2f"%(seq, self.oneletter[offset:offset+len(seq)],
                                                                   total_f, total_r,  maximum,
                                                                   max([total_f,total_r])/self.maxscore)
            if total_f > threshold and total_f > total_r:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_f)
                matches.append(seq[offset:offset+self.width])
            elif total_r > threshold:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_r)
                matches.append(seq[offset:offset+self.width])
        return matches,endpoints,scores                

    def mask_seq(self,seq):
        """return a copy of input sequence in which any regions matching m are
        replaced with strings of N's """
        masked = ''
        matches, endpoints, scores = self.scan(seq)
        cursor = 0
        for start, stop in endpoints:
            masked = masked + seq[cursor:start] + 'N'*self.width
            cursor = stop+1
        masked = masked + seq[cursor:]
        return masked

    def masked_neighborhoods(self,seq,flanksize):
        """chop up the input sequence into regions surrounding matches to m.
        Replace the subsequences that match the motif with N's."""
        ns = self.seq_neighborhoods(seq,flanksize)
        return [self.mask_seq(n) for n in ns]

    def seq_neighborhoods(self,seq,flanksize):
        """chop up the input sequence into regions surrounding matches to the
        motif."""
        subseqs = []
        matches, endpoints, scores = self.scan(seq)
        laststart, laststop = -1, -1
        for start, stop in endpoints:
            curstart, curstop = max(0,start-flanksize), min(stop+flanksize,len(seq))
            if curstart > laststop:
                if laststop != -1:
                    subseqs.append(seq[laststart:laststop])
                laststart, laststop = curstart, curstop
            else:
                laststop = curstop
        if endpoints: subseqs.append(seq[laststart:laststop])
        return subseqs

    def __sub__(self,other):
        pass
        """Overloads the '-' operator to compute the Euclidean distance between
        probability matrices motifs of equal width."""
        if type(other) != type(self):
            print "computing distance of unlike pssms (types %s, %s)"%(
                type(other),type(self))
            print 'First: %s'%other
            print 'Self:  %s'%self
            sys.exit(1)
        if other.width != self.width:
            print "computing distance of unlike pssms (width %d != %d)"%(
                other.width,self.width)
            sys.exit(1)
        D = 0
        FABS = math.fabs
        POW  = math.pow
        for L in self.logP[0].keys():
            for i in range(self.width):
                D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
                #D = D + FABS( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]))
                #D = D + FABS(self.logP[i][L] - other.logP[i][L])
        return math.sqrt(D)

    def maskdiff(self,other):
        """a different kind of motif comparison metric.  See THEME paper for
        details"""
        return maskdiff(self,other)

    def maxdiff(self):
        """compute maximum possible Euclidean distance to another motif.  (For
        normalizing?)"""
        POW  = math.pow
        D = 0
        for i in range(self.width):
            _min = 100
            _max = -100
            for L in ACGT:
                val = POW(2,self.logP[i][L])
                if   val > _max:
                    _max  = val
                    _maxL = L
                elif val < _min:
                    _min  = val
                    _minL = L
            for L in ACGT:
                if L == _minL:
                    delta = 1-POW(2,self.logP[i][L])           #1-val
                    D = D + delta*delta
                else:
                    D = D + POW( POW(2,self.logP[i][L]), 2)    #0-val
        return math.sqrt(D)
                
    def revcomp(self):
        """return reverse complement of motif"""
        return revcompmotif(self)
    def trimmed(self,thresh=0.1):
        """return motif with low-information flanks removed.  'thresh' is in bits."""
        for start in range(0,self.width-1):
            if self.bits[start]>=thresh: break
        for stop  in range(self.width,1,-1):
            if self.bits[stop-1]>=thresh: break
        m = self[start,stop]
        return m
    def bestseqs(self,thresh=None):
        """return all k-mers that match motif with a score >= thresh"""
        if not thresh:
            if self._bestseqs:
                return self._bestseqs
        if not thresh: thresh = 0.8 * self.maxscore
        self._bestseqs = bestseqs(self,thresh)
        return self._bestseqs
    def emit(self,prob_min=0.0,prob_max=1.0):
        """consider motif as a generative model, and have it emit a sequence"""
        if not self.cumP:
            for logcol in self.logP:
                tups = []
                for L in ACGT:
                    p = math.pow(2,logcol[L])
                    tups.append((p,L))
                tups.sort()
                cumu = []
                tot  = 0
                for p,L in tups:
                    tot = tot + p
                    cumu.append((tot,L))
                self.cumP.append(cumu)
        s = []
        #u = random()+0.01 #Can make higher for more consistent motifs
        for cumu in self.cumP:
            u = (prob_max-prob_min)*random() + prob_min
            #u = random()+0.01 #Can make higher for more consistent motifs
            last = 0
            for p,L in cumu:
                if last < u and u <= p:
                    letter = L
                    break
                else: last = p
#           print L,'%8.4f'%u,cumu
            s.append(L)
        #print ''.join(s)
        return ''.join(s)
            
                
    def random_kmer(self):
        """generate one of the many k-mers that matches the motif.  See m.emit()
        for a more probabilistic generator"""
        if not self._bestseqs: self._bestseqs = self.bestseqs()
        seqs   = self._bestseqs
        pos = int(random() * len(seqs))
        print 'Random: ',self.oneletter,seqs[pos][1]
        return seqs[pos][1]

    def __getitem__(self,tup):
        pass
        """
        m.__getitem__(tup) -- Overload m[a,b] to submotif.  Less pythonish than [:], but more reliable
        """
        if len(tup) != 2:
            print "Motif[i,j] requires two arguments, not ",tup
        else:
            beg, end = tup[0], tup[1]
            return submotif(self,beg,end)
    def __getslice__(self,beg,end):
        pass
        """
        m.__getslice__(,beg,end) -- Overload m[a:b] to submotif.
        """
        if beg >= end:
            #Probably python converted negative idx.  Undo
            beg = beg - self.width
        return submotif(self,beg,end)
    def __add__(self,other):
        pass
        """
        m.__add__(other) -- Overload  '+' for concatenating motifs
        """
        return merge(self,other,0)
    def __len__(self):
        pass
        """
        m.__len__()  -- Overload len(m) to return width
        """
        return self.width
    def shuffledP(self):
        """
        m.shuffledP() -- Generate motif in which probability matrix has been shuffled.
        """
        return shuffledP(self)
    def copy(self):
        """return a 'deep' copy of the motif"""
        a = Motif()
        a.__dict__ = self.__dict__.copy()
        return a

    def random_diff_avestd(self,iters=5000):
        """see modules' random_diff_avestd"""
        return random_diff_avestd(self,iters)
    def bogus_kmers(self,count=200):
        """Generate a faked multiple sequence alignment that will reproduce the
        probability matrix."""

        POW  = math.pow
        #Build p-value inspired matrix
        #Make totals cummulative:
        # A: 0.1 C: 0.4 T:0.2 G:0.3
        #                            ->  A:0.0 C:0.1 T:0.5 G:0.7  0.0
        
        #Take bg into account:
        # We want to pick P' for each letter such that:
        #     P'/0.25  = P/Q
        # so  P'       = 0.25*P/Q
        
        m = []
        for i in range(self.width):
            _col = []
            tot   = 0.0
            for L in ACGT:
                _col.append( tot )
                tot = tot + POW(2,self.logP[i][L]) * 0.25 / self.background[L]
            _col.append(tot)
            #Renormalize
            for idx in range(len(_col)):
                _col[idx] = _col[idx] / _col[-1]
            m.append(_col)

        for p in range(0): #Was 5
            for i in range(len(m)):
                print '%6.4f  '%m[i][p],
            print

        seqs=[]
        for seqnum in range(count+1):
            f = float(seqnum)/(count+1)
            s = []
            for i in range(self.width):
                for j in range(4):
                    if (m[i][j] <= f and f < m[i][j+1]):
                        s.append(ACGT[j])
                        break
            seqs.append(''.join(s))

        del(seqs[0])
        #for i in range(count):
        #    print ">%3d\n%s"%(i,seqs[i])

        return seqs


def minwindowdiff(M1,M2,overlap=5,diffmethod='diff'):
    #Alternate method: maskdiff, infomaskdiff
    if type(M1) != type(M2):
        print "Error: Attempted to compute alignment of objects that are not both Motifs"
        print "       types %s: %s  and %s: %s"%(M1,type(M1),M2,type(M2))
        sys.exit(1)

    if M1.width <= M2.width: A = M1; Borig = M2
    else:                    A = M2; Borig = M1
    wA = A.width
    wB = Borig.width
    O  = overlap

    if   diffmethod == 'diff':
        diff_fcn = diff
    elif diffmethod == 'maskdiff':
        diff_fcn = maskdiff
    elif diffmethod == 'infomaskdiff':
        diff_fcn = infomaskdiff
        
    mindiff = 1000
    #print 'minwindodebug    wA ', wA, 'wB ', wB, 'O ', O, 'wA-0', wA-O, 'wB-O', wB-O
    for Astart in range(wA-O+1):
        subA = A[Astart:Astart+O]
        for B in [Borig, Borig.revcomp()]:
            for Bstart in range(wB-O+1):
                subB = B[Bstart:Bstart+O]
                mindiff = min(mindiff, diff_fcn(subA,subB))
                #print 'minwindodebug     ',subA, subB, diff_fcn(subA,subB)
    return mindiff
    

def minaligndiff(M1,M2,overlap=5,diffmethod='diff'):
    #Alternate method: maskdiff, infomaskdiff
    if type(M1) != type(M2):
        print "Error: Attempted to compute alignment of objects that are not both Motifs"
        print "       types %s: %s  and %s: %s"%(M1,type(M1),M2,type(M2))
        sys.exit(1)

    if M1.width <= M2.width:
        A = M1; Borig = M2
        switch = 0
    else:
        A = M2; Borig = M1
        switch = 1
    wA = A.width
    wB = Borig.width
    O  = overlap

    '''
    Here is the figure to imagine:
       012345678901234567890   wA: 6  Bstart: 6-3     = 3
         A         (A)         wB: 11 Bstop:  6+11-3-1= 13
       ------     %%%%%%        O: 3  lastA:  6+11-3-3= 11
          -----------
          |O|  B
    '''

    if   diffmethod == 'diff':
        diff_fcn = diff
    elif diffmethod == 'maskdiff':
        diff_fcn = maskdiff
    elif diffmethod == 'infomaskdiff':
        diff_fcn = infomaskdiff
    
    Bstart = wA-O
    Bstop  = wA+wB-O-1
    lastA  = wA+wB-O-O
    Dmin = 1000
    Dmins=[]
    #print A
    #print '%s%s'%(' '*Bstart,Borig)
    for B in [Borig, Borig.revcomp()]:
        for start in range(0,lastA+1):
            Bpos = []
            Apos = []
            for offset in range(wA):
                abs = start+offset
                if abs >= Bstart and abs <= Bstop:
                    Apos.append(offset)
                    Bpos.append(abs-Bstart)
            subA = A[min(Apos),max(Apos)+1]
            subB = B[min(Bpos),max(Bpos)+1]
            #print '%s%s\n%s%s  %f'%(
            #    ' '*start, subA,
            #    ' '*start, subB,   diff_fcn(subA,subB))
            if switch: _diff = diff_fcn(subB,subA)
            else:      _diff = diff_fcn(subA,subB)
            Dmin = min(Dmin, _diff)
    return Dmin
    
'''
To compare 2 motifs of the same width, there are these five functions:

m1 - m2            - Euclidean Distance (sqrt(sum_col(sum_row)))
diff(m1,m2)        - psuedo-Euclidean (sum_col(sqrt(norm(sum_row)))/#col
maskdiff(m1,m2)    - diff, but excluding positions with "N" in m2
infomaskdiff(m1,m2)- diff, but scaling distance by normalized
     information content at each position in m2.
diverge(m1,m2)     - Mutual information sum[p log (p/q)]

**Note that maskdiff, infomaskdiff, and diverge are not symmetric functions

To compare 2 motifs of different widths, there is the function:

minaligndiff(M1,M2,overlap=5,diffmethod='diff')

this does a 'sliding' comparison of two motifs and reports the minimum
distance over all alignments.  overlap refers to the minumum overlap
required while sliding.  Below, overlap is '2'.  The default is '5'.

      ------
          -----------

You can optionally specify the distance metric as a text string.
The default is 'diff'.

'''


def diff(self,other):
    """psuedo-Euclidean (sum_col(sqrt(norm(sum_row)))/#col"""
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)
    POW     = math.pow
    Dtot    = 0
    for i in range(self.width):
        '''Computes distance'''
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        Dtot = Dtot + math.sqrt(D)/math.sqrt(2.0)
    return Dtot/self.width
    

def maskdiff(self,other):
    """diff, but excluding positions with 'N' in m2. Return pseudo-Euclidean
    distance, but only include columns that are not background."""
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    Dtot = 0
    POW  = math.pow
    NEAR0= lambda x:(-0.01 < x and x < 0.01)
    divisor = 0
    for i in range(self.width):
        nearcount = 0

        '''Implements mask'''
        for L in ACGT:
            diff = POW(2,other.logP[i][L]) - other.background[L]
            if NEAR0(diff): nearcount = nearcount + 1
        if nearcount == 4:
            #print 'Skipping position %d :'%i,other.logP[i]
            continue

        '''Computes distance'''
        divisor = divisor + 1
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        Dtot = Dtot + math.sqrt(D)/math.sqrt(2.0)
    return Dtot/divisor

def infomaskdiff(self,other):
    """Return pseudo-Euclidean distance, but scale column distance by
    information content of "other".  Used by THEME"""
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    maxbits = math.log( 1.0/min(other.background.values()) ) / math.log(2.0)
    '''or... alternatively'''
    #print maxbits, max(other.bits)
    #print other.bits
    maxbits = max(other.bits)
    if maxbits < 0.1:  #'''There is nothing important here'''
        return 1
    
    Dtot    = 0
    POW     = math.pow
    divisor = 0
    '''Computes distance'''
    for i in range(self.width):
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        col_dist  = math.sqrt(D)/math.sqrt(2.0)
        col_scale = other.bits[i]/maxbits
        divisor = divisor + col_scale
        Dtot = Dtot + col_dist*col_scale
    return Dtot/divisor

def diverge(self,other):
    """Yet another distance metric"""
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    Dtot = 0
    POW  = math.pow
    LOG2 = lambda x:math.log(x)/math.log(2.0)
    NEAR0= lambda x:(-0.01 < x and x < 0.01)
    divisor = 0
    for i in range(self.width):
        nearcount = 0

        '''Implements mask'''
        for L in ACGT:
            diff = POW(2,other.logP[i][L]) - self.background[L]
            if NEAR0(diff): nearcount = nearcount + 1
        if nearcount == 4:
            #print 'Skipping position %d :'%i,other.logP[i]
            continue

        '''Computes distance'''
        divisor = divisor + 1
        D = 0
        for L in ACGT:
            Pself = POW(2, self.logP[i][L])
            Pother= POW(2,other.logP[i][L])
            D = D + Pself * LOG2(Pself/Pother)
        Dtot = Dtot + D
    return Dtot/divisor



def bestseqs(motif,thresh, seq='',score=0,depth=0,bestcomplete=None,SEQS=[]):
    """This function returns a list of all sequences that a motif could
    match match with a sum(log-odds) score greater than thresh."""
    if depth == 0:
        SEQS = []  #Must be a python 2.1 bug. I shouldn't have to do this
    if not bestcomplete:
        M = motif
        maxs = []
        for i in range(M.width):
            bestj = 'A'
            for j in ['C', 'G', 'T']:
                if M.ll[i][j] > M.ll[i][bestj]:
                    bestj = j
            maxs.append(M.ll[i][bestj])
        bestcomplete = []
        for i in range(M.width):
            tot = 0
            for j in range(i,M.width):
                tot = tot + maxs[j]
            bestcomplete.append(tot)
    if depth == motif.width:
        if score > thresh:
            SEQS.append((score,seq))
        #if len(SEQS) > 2000:
        #    thresh = 1000.0 # Return Early, You don't really want all these sequences, do you?
        return
    if depth==-1:
        print '# %-10s %6.3f %6.3f %2d'%(seq, score, bestcomplete[depth], depth)
    if score + bestcomplete[depth] < thresh: return
    #if depth > 0 and len(SEQS) > 2000:
    #    return
    for L in ACGT:
        newseq   = seq + L
        newscore = score + motif.ll[depth][L]
        bestseqs(motif,thresh,newseq,newscore,depth+1,bestcomplete,SEQS)
    if depth == 0:
        SEQS.sort()
        SEQS.reverse()
        return SEQS

def seqs2fasta(seqs,fasta_file = ''):
    """
    seqs2fasta(seqs,fasta_file = '') -- Dumps a Fasta formatted file of sequences,
    keyed by the sequence itself::

      >ACTTTTTGTCCCA
      ACTTTTTGTCCCA
      >ACTTTTGGGGCCA
      ACTTTTGGGGCCA
        ...

    """
    if not fasta_file:
        fasta_file = tempfile.mktemp()
    FH = open(fasta_file,'w')
    for i in range(len(seqs)):
        FH.write(">%d\n%s\n"%(i,seqs[i]))
    FH.close()
    return fasta_file

def top_nmers(N,seqs,with_counts = 0,purge_Ns = ''):
    """Assemble list of all nmers (kmers) with width 'N' from supplied sequences.
    Option with_counts returns list of (kmer, count) tuples instead.  Purge N's
    ignores kmers containing N's.  """
    Nmers = {}
    revcompTBL = string.maketrans("AGCTagctnN", "TCGAtcganN")
    for seq in seqs:
        for i in range(len(seq)-N+1):
            Nmer = seq[i:i+N]
            if purge_Ns:
                if Nmer.find('N') >= 0: continue
            _t = list(Nmer.translate(revcompTBL))
            _t.reverse()
            NmerRC = ''.join(_t)   # _t used until here to revese comp seq
            _t = [Nmer, NmerRC]
            _t.sort()
            NmerKey = _t[0]        # _t used until here to get alphabetically first seq
            if Nmers.has_key(NmerKey):
                Nmers[NmerKey] = Nmers[NmerKey] + 1
            else:
                Nmers[NmerKey] = 1
    sorted = Nmers.keys()
    sorted.sort(lambda x,y,D=Nmers:cmp(D[y],D[x]) or cmp(x,y))
    #for i in range(10):
    #    print "# %2d  %s %d"%(i,sorted[i],Nmers[sorted[i]])
    if with_counts:
        return zip(sorted,map(lambda x,N=Nmers:N[x], sorted))
    else:
        return sorted

def m_matches(seqs,wmer,m):
    """Returns list of all kmers among sequences that have at most
    m mismatches to the supplied wmer (kmer)."""
    matches = []
    width = len(wmer)
    for (nmer, count) in top_nmers(width,seqs,'with counts'):
        match = 0
        for i in range(width):
            if nmer[i] == wmer[i]:
                match = match+1
        if match >= m:
            for i in range(count):
                matches.append(nmer)
    return matches

def compare_seqs(s1, s2):
    pass
    """
    compare_seqs(s1, s2) 
    """
    if len(s1) > len(s2):
        long  = s1
        short = s2
    else:
        long  = s2
        short = s1
    (maxcount,max_i) = (0,0)
    for i in range(len(long)-len(short)+1):
        idcount_f = 0
        idcount_r = 0
        for j in range(len(short)):
            if short[j] == long[i+j]:
                idcount_f = idcount_f + 1
            if short[-(j+1)] == revcomp[long[i+j]]:
                idcount_r = idcount_r + 1
        if (idcount_f > maxcount and idcount_f >= idcount_r):
            maxcount = idcount_f
            max_i    = i
        elif (idcount_r > maxcount):
            maxcount = idcount_r
            max_i    = i
        #print i,j,idcount_f,idcount_r,maxcount
    maxfrac = float(maxcount) / len(short)
    print maxfrac,maxcount,len(short)
    return maxfrac,short,long[max_i:max_i+len(short)]

def shuffle_bases(m):
    """return a new motif object in which the probabilities are randomly
    re-assigned to different letters at the same position."""
    C = []
    letts = list('ACGT')
    for i in range(m.width):
        D = {}
        vals = m.counts[i].values()
        shuffle(vals)
        for i in range(4):
            D[letts[i]] = vals[i]
        C.append(D)
    n = Motif()
    #n.__dict__ = m.__dict__.copy() #May copy too much information (cached diff information, etc...)
    n.compute_from_counts(C)
    return n

def random_diff_avestd(motif,iters=5000):
    """Return the average & stddev distance ('diff') between a
    motif and "iters" random motifs of the same width."""
    w = motif.width
    vals = []
    for i in range(iters):
        vals.append(motif - Random_motif(w))
    return avestd(vals)

def random_motif(w):
    """Generate a random motif of width w.  Each position will have a dominant
    letter with probability around 0.91."""
    C = []
    for i in range(w):
        D = {}
        tot = 0
        p = int(random.random() * 4)
        Lup = ACGT[p]
        for L in ACGT:
            D[L] = 0.1
            tot = tot + 0.001
        D[Lup] = D[Lup] + 1
        for L in ACGT:
            D[L] = D[L]/tot
        C.append(D)
    m = Motif()
    m.compute_from_counts(C)
    return m

def toDict(M):
    pass
    '''
    toDict(M) -- Convert a 2D array to a list of dictionaries (which is how the motif object
                 stores information internally).  Assumes M entries are in alphabetical order (ACGT)
    '''
    if type(M[0]) == type(0.0):
        return toDictVect(M)
    else:
        a = []
        for i in range(len(M)):
            a.append(toDictVect(M[i]))
        return a
        
def toDictVect(V):
    pass
    """
    toDictVect(V) -- Convert a 1D vector to a dictionary of DNA letters.  Assumes values
    in V are in alphabetical order (ACGT).
    """
    D = {}
    for L,i in (('A',0), ('C',1), ('G',2), ('T',3)):
        D[L]=V[i]
    return D

def submotif(self,beg,end):
    """**Deprecated** Use slice functionality (m[2:4]) instead.
    
    Utility function
    for extracting sub-motifs and padding motifs."""
    bg = self.background.copy()
    P = []

    #Determine if any 'zeros' should be added at begining
    #because the user has specified a negative beg index
    for i in range(beg,0):
        P.append(bg.copy())

    #Copy relevant content of motif
    start = max(beg,0)
    stop  = min(end,self.width)
    for i in range(start,stop):
        D = {}
        for L in ACGT:
            D[L] = math.pow(2.,self.logP[i][L])
        P.append(D)

    #Determine if any 'zeros' should be added at the end
    #because the user has specified a width too large
    for i in range(self.width,end):
        P.append(bg.copy())

    #print "BEG, END", beg,end
    #for i in range(beg,end):
    #    print i,P[i]

    #Build the Motif
    M = copy.deepcopy(self)
    #M = Motif(None,bg.copy())
    M.compute_from_counts(P)
    M.source = self.source
    return M
                
def shuffledP(self):
    """Construct a motif in which the letter distributions are preserved but
    are reassigned to rondom positions in the motif."""
    bg = self.background.copy()
    P = []

    #Copy relevant content of motif
    for i in range(0,self.width):
        D = {}
        _s = ACGT[:]
        shuffle(_s)
        for L,_L in zip(ACGT,_s):
            D[L] = math.pow(2.,self.logP[i][_L])
        P.append(D)

    #Build the Motif
    M = copy.deepcopy(self)
    #M = Motif(None,bg.copy())
    M.compute_from_counts(P)
    M.source = self.source
    return M

def revcompmotif(self):
    """Construct the reverse complement of the motif.  Use m.revcomp() member
    function instead."""
    bg = self.background.copy()
    P = []

    for i in range(self.width):
        D = {}
        for L in ACGT:
            D[L] = math.pow(2.,self.logP[self.width-i-1][revcomp[L]])
        P.append(D)

    #Build the Motif
    M = copy.deepcopy(self)
    M.compute_from_counts(P)
    return M
        

def sum(motifs,weights=[]):
    """Perhaps better called 'average'.  Constructs a motif by averaging the
    probabilities at each position of the (pre-aligned) input motifs.  Optional
    weights can be assigned, and must be in the same order as the motifs. 
    """
    if not weights:
        weights = [1.0] * len(motifs)
    tot = 0.0
    for w in weights: tot=tot+float(w)
    weights = [(w/tot) for w in weights]
    C = []
    for c in motifs[0].fracs:
        D = {}
        for L in ACGT: D[L] = 0.0
        C.append(D)
    for m,w in zip(motifs,weights):
        for i in range(m.width):
            for L in ACGT:
                C[i][L] = C[i][L] + m.fracs[i][L]*w
    motif = Motif_from_counts(C,0.0,bg=motifs[0].background)
    return motif.trimmed()


def giflogo(motif,id,title=None,scale=0.8):
    """Interface to the 'weblogo/seqlogo' perl
    scripts that generate colorful sequence logos
    """
    return seqlogo(motif,id,title,scale,format='GIF')


seqlogo_formats = ('GIF','PDF','EPS','PNG')
illegal_fn_chars = '&;/ ()'
fn_trans = string.maketrans(illegal_fn_chars,'_'*len(illegal_fn_chars))
def seqlogo(motif,motif_id,title=None,scale=0.8,img_format='GIF') :
    """Interface to the'weblogo/seqlogo' perl scripts that generate colorful
    sequence logos.  Available formats are %s.  Replaces illegal filename
    characters in *id* parameter (i.e. '%s') with underscores when writing
    to file.  The executable *seqlogo* must be on your path.
    """%(seqlogo_formats,illegal_fn_chars)
    #SEQLOGO = TAMOpaths.weblogodir + 'seqlogo'
    #TAMOpaths.CHECK(SEQLOGO,'','Weblogo/Seqlogo')
    kmers   = motif.bogus_kmers(100)
    width   = float(len(kmers[0]) )
    height  = float(4)
    m       = motif
    width, height = width*scale, height*scale
    tmp     = tempfile.mktemp() + '.fsa'
    if title is None:
        title = motif_id

    if img_format.upper() not in seqlogo_formats :
        raise MotifToolsException('seqlogo requires one of %s'%seqlogo_formats)

    seqs2fasta(kmers,tmp)
    fn = id.translate(fn_trans)
    cmd = 'seqlogo -F %s -acpY -w%d -h%d -k 1 -M -f %s -o %s -t "%s" '%(
          img_format.upper(), width, height, tmp, fn, title)

    call(cmd,shell=True)
    return "%s.%s"%(fn,img_format.lower())


def merge(A,B,overlap=0):
    """**Deprecated** Use the '+' operator instead.
    
    Used for concatenating motifs into a new motif, allowing for the averaging
    of overlapping bases between them.
    """
    if (overlap < 0 or overlap > A.width or overlap >B.width):
        print 'Cannot overlap %s with %s by %d bases'%(A.oneletter,B.oneletter,overlap)
        return None

    #Build Probability matrix.  Width will be A.width + B.width - overlap
    w = A.width + B.width - overlap

    P = []
    #Make a copy of A's probabilities into P
    for i in range(A.width):
        D = {}
        logP = A.logP[i]
        for L in logP.keys():
            D[L] = math.pow(2,logP[L])
        P.append(D)
    #Add B's first 'overlap' probabilities to last 'overlap' probabilities of P
    for i in range(overlap):
        logP = B.logP[i]
        Pidx = len(P)-overlap+i
        _tot = 0
        for L in logP.keys():
            P[Pidx][L] = (P[Pidx][L] + math.pow(2,logP[L])) / 2.0
            P[Pidx][L] = max(P[Pidx][L],math.pow(2,logP[L]))
            _tot = _tot + P[Pidx][L]
        for L in logP.keys():
            P[Pidx][L] = P[Pidx][L] / _tot
    #Append B's remaining probabilites to P
    for i in range(overlap,B.width):
        D = {}
        logP = B.logP[i]
        for L in logP.keys():
            D[L] = math.pow(2,logP[L])
        P.append(D)
        
    #Build a motif
    M = Motif(None,A.background.copy())
    M.source = A.source,B.source
    M.compute_from_counts(P)
    return M

def avestd(vals):
    """return an (average, stddev) tuple computed from the supplied list of values"""
    (sum, sum2) = (0.,0.)
    N = float(len(vals))
    for val in vals:
        sum  = sum  + float(val)
        sum2 = sum2 + float(val)*float(val)
    if N == 1:
        ave = sum
        std = 0
    else:
        ave = sum /  N
        std = math.sqrt( (sum2-(N*ave*ave)) / (N-1.0) )
    return ave,std


def load(filename):
    """load a 'TAMO'-formatted motif file"""
    FID = open(filename,'r')
    lines = FID.readlines()
    FID.close()
    motifs   = []
    seedD    = {}
    seedfile = ''
    for i in range(len(lines)):
        if lines[i][0:10] == 'Log-odds matrix'[0:10]:
            w = len(lines[i+1].split())-1
            ll = []
            for pos in range(w):
                ll.append({})
            for j in range(0,4):
                toks = lines[i+j+2].split()
                L = toks[0][1]
                for pos in range(w):
                    ll[pos][L] = float(toks[pos+1])
            m = Motif_from_ll(ll)
            motifs.append(m)
        if lines[i][0:6] == 'Motif '[0:6]:
            toks =  lines[i].split()
            motifs[-1].nseqs    = float(re.sub('[\(\)]','',toks[3]))
            motifs[-1].totalbits= float(toks[5])
            motifs[-1].MAP      = float(toks[7])
            motifs[-1].seeddist = float(toks[9])
            motifs[-1].seednum  = int(toks[10][0:-1])
            motifs[-1].pvalue   = math.pow(10,-float(toks[12]))

            if 'ch:' in toks:
                _idx = toks.index('ch:')
                motifs[-1].church = math.pow(10,-float(toks[_idx+1]))
            if 'Es:' in toks:
                _idx = toks.index('Es:')
                motifs[-1].E_site = math.pow(10,-float(toks[_idx+1]))
            if 'x2:' in toks:
                _idx = toks.index('x2:')
                motifs[-1].E_chi2 = math.pow(10,-float(toks[_idx+1]))
            if 'Eq:' in toks:
                _idx = toks.index('Eq:')
                motifs[-1].E_seq = math.pow(10,-float(toks[_idx+1]))
            if 'mn:' in toks:
                _idx = toks.index('mn:')
                motifs[-1].MNCP = float(toks[_idx+1])
            if 'f:' in toks:
                _idx = toks.index('f:')
                motifs[-1].frac = float(toks[_idx+1])
            if 'Ra:' in toks:
                _idx = toks.index('Ra:')
                motifs[-1].ROC_auc = float(toks[_idx+1])
            if 'cR:' in toks:
                _idx = toks.index('cR:')
                motifs[-1].CRA     = float(toks[_idx+1])
            if 'Cf:' in toks:
                _idx = toks.index('Cf:')
                motifs[-1].Cfrac   = float(toks[_idx+1])
            if 'k:' in toks:
                _idx = toks.index('k:')
                motifs[-1].kellis  = float(toks[_idx+1])

            if 'b:' in toks:
                _idx = toks.index('b:')
                motifs[-1].numbound = int(toks[_idx+1])
            if 'nG:' in toks:
                _idx = toks.index('nG:')
                motifs[-1].nummotif = int(toks[_idx+1])
            if 'bn:' in toks:
                _idx = toks.index('bn:')
                motifs[-1].numboundmotif = int(toks[_idx+1])



        if lines[i][0:10] == 'Threshold: '[0:10]:
            toks =  lines[i].split()
            motifs[-1].threshold= float(toks[1])
        if lines[i][0:5] == 'Seed '[0:5]:
            toks = lines[i].split()
            id = int(toks[1][0:-1])  #'10:' -> '10'
            seedD[id] = toks[2]
        if lines[i][0:7] == 'Source: '[0:7]:
            motifs[-1].source = lines[i][7:].strip()
        if lines[i][0:6] == 'Gamma: '[0:6]:
            motifs[-1].gamma = float(lines[i][6:])
        if lines[i][0:6] == 'Evalue: '[0:6]:
            motifs[-1].evalue = float(lines[i][7:].strip())
        if lines[i][0:22]=='Program specific score: '[0:22]:
            tempprogscore=lines[i][23:].split(":");

            for i in range(len(tempprogscore)):
                tempprogscore[i]=tempprogscore[i].strip()

            if len(tempprogscore)>1:
                try:
                    tempprogscore[1]=float(tempprogscore[1])
                except ValueError:
                    tempprogscore[1]=tempprogscore[1]
                motifs[-1].progscore=tempprogscore

        if lines[i][0:10] == 'fasta file:'[0:10]:
            parts=lines[i].strip().split()
            motifs[-1].dataset, motifs[-1].beta, motifs[-1].bgfile = \
                        parts[2],float(parts[4]), parts[7]

    if lines[i][0:21]=='classification error: '[0:21]:
        motifs[-1].cverror=float(lines[i][22:].strip())
    if lines[i][0:20]=='SVM match threshold: '[0:20]:
        motifs[-1].match_thresh=float(lines[i][21:].strip())
    if lines[i].find('Using')>=0 and lines[i].find('as seeds')>=0:
            '''#Using all (132) motifs in SLT_081503.seeds as seeds:'''
            seedfile = lines[i].split()[-3]
    for i in range(len(motifs)):
        if seedfile: motifs[i].seedfile = seedfile
        seednum = motifs[i].seednum
        if seedD.has_key(seednum):
            motifs[i].seedtxt = seedD[seednum]
    return motifs
    
def save_motifs(motifs,filename,kmer_count=20):
    """Save list of motifs as a 'TAMO'-formatted motif file to the specificied file.
    optional kmer_count specificies how many sequences to include in the printed
    multiple sequence alignment that recapitulates the probability matrix."""
    try :
        print_motifs(motifs,kmer_count,f=filename)
    except:
        print '!-- Error saving motifs to %s'%filename
        raise
    
def print_motif(motif,kmer_count=20,istart=0,f=None):
    """Print a motif in the 'TAMO'-format.  istart specificies the motif number, and 
    optional kmer_count specificies how many sequences to include in the printed
    multiple sequence alignment that recapitulates the probability matrix. """
    print_motifs([motif],kmer_count,istart)
    sys.stdout.flush()

def print_motifs(motifs,kmer_count=20,istart=0,f=None):
    """Print list of motifs as a 'TAMO'-formatted motif file to the specificied file.
    Optional kmer_count specificies how many sequences to include in the printed
    multiple sequence alignment that recapitulates the probability matrix.
    istart specifies number from which to begin motif ids."""

    # handle f input cases
    if f is None :
        f = sys.stdout
    elif isinstance(f,str) :
        f = open(f,'w')

    i = istart-1
    for m in motifs:
        i = i + 1
        print >>f,  "Log-odds matrix for Motif %3d %s"%(i,m)
        m._print >>f, _ll()
        #print >>f,  "Probability matrix for Motif %3d %s"%(i,m)
        #m._print >>f, _p()
        print >>f,  "Sequence Logo"
        m._print >>f, _bits()
        for newprop in ('gamma', 'church', 'E_site', 'E_seq', 'E_chi2', 'realpvalue',
                        'kellis', 'MNCP', 'ROC_auc', 'CRA', 'Cfrac', 'frac', 'binomial'):
            if not m.__dict__.has_key(newprop):   #Kludge to deal w/ old shelves
                m.__dict__[newprop] = None
        if m.seedtxt:  print >>f,  "Seed: %3d %s"%(i,m.seedtxt)
        if m.gamma:    print >>f,  "Gamma: %7.5f"%m.gamma
        if m.evalue != None: print >>f,  'Evalue: %6.3e'%m.evalue
        if m.progscore is not None :
            printableProgscore=(m.progscore[0],str(m.progscore[1]))
            print >>f,  'Program specific score: '+ ": ".join(printableProgscore)

        if m.family:   print >>f,  "Family: ",m.family
        if m.source:   print >>f,  "Source: ",m.source
        if m.dataset:  print >>f,  "fasta file: %s beta: %f background sequences: %s"%(m.dataset,m.beta,m.bgfile)
        if m.match_thresh: print >>f,  "SVM match threshold: ",m.match_thresh
        if m.cverror:  print >>f,  "classification error: ",m.cverror
        #Motif   0 NGAGGGGGNN (0)            (Bits:   8.24   MAP:   6.53   D:  0.21  0)  Enr: 54.000 
        print >>f,  "Motif %3d %-25s (Bits: %5.2f  MAP: %5.2f   D: %5.3f  %2d) E: %6.3f"%(
            i, m, m.totalbits, m.MAP, m.seeddist, m.seednum, nlog10(m.pvalue)),
        if m.binomial!=None:  print >>f,  ' Bi: %5.2f'%nlog10(m.binomial),
        if m.church != None:  print >>f,  ' ch: %5.2f'%nlog10(m.church),
        if m.frac   != None:  print >>f,  ' f: %5.2f'%(m.frac),
        if m.E_site != None:  print >>f,  ' Es: %5.2f'%nlog10(m.E_site),
        if m.E_seq != None:  print >>f,  ' Eq: %5.2f'%(nlog10(m.E_seq)),
        if m.MNCP   != None:  print >>f,  ' mn: %5.2f'%(m.MNCP),
        if m.ROC_auc!= None:  print >>f,  ' Ra: %6.4f'%(m.ROC_auc),
        if m.E_chi2 != None:
            if m.E_chi2 == 0: m.E_chi2=1e-20
            print >>f,  ' x2: %5.2f'%(nlog10(m.E_chi2)),
        if m.CRA    != None:  print >>f,  ' cR: %6.4f'%(m.CRA),
        if m.Cfrac  != None:  print >>f,  ' Cf: %6.4f'%(m.Cfrac),
        if m.realpvalue != None: print >>f,  ' P: %6.4e'%(m.realpvalue)
        if m.kellis != None:  print >>f,  ' k: %5.2f'%(m.kellis),
        try:
            if m.numbound      :  print >>f,  ' b: %3d'%(m.numbound),
            if m.nummotif      :  print >>f,  ' nG: %3d'%(m.nummotif),
            if m.numboundmotif :  print >>f,  ' bn: %3d'%(m.numboundmotif),
        except: pass
        print >>f, ''

        _max = m.maxscore
        m.maxscore = -100
        if kmer_count >= 0:
            seqs = m.bogus_kmers(kmer_count)
        else:
            seqs = m.seqs

        for seq in seqs:
            print >>f,  seq,i,m.scan(seq)[2][0]

        m.maxscore = _max
        print >>f,  '*'*m.width
        print >>f,  "MAP Score: %f"%(m.MAP)

def nlog10(x,min=1e-323):
    """returns -log10(x) with a maximum default value of 323."""
    if x < min: x=min
    try:
        return math.fabs(math.log(x)/math.log(10))
    except:
        return 0

def txt2motifs(txt,VERBOSE=1):
    """Convert a text string into a list of motifs:
    Examples:

    'TGASTCA,GAATC'      --> 2 motifs from ambiguity codes
    'results.tamo'       --> All motifs in TAMO-format file
    'results.tamo:34,45' --> Motifs 34 and 45 in TAMO-format file
    'results.pickle'     --> All motifs in pickle (list or dict of Motifs)
    'results.pickle%GAL4 --> 'GAL4' entry in results.pickle dictionary
    'results.pickle:34,45 -> Motifs 34 and 45 in results.pickle list
    """
    motifs = []
    exists = os.path.exists
    toks   = txt.split(':')
    if exists(toks[0]):               #It's a file!!
        fname = toks[0]
        if fname.find('.pickle') > 0: #It's a pickle!!
            return pickletxt2motifs(toks)
        else:                         #It's a "Motif" file!!
            if VERBOSE:
                print "# Loading motif from %s"%fname
            allmotifs = load(fname)
        if len(toks) == 1: motifs = allmotifs
        else:
            idxs   = [int(x) for x in toks[1].split(',')]
            motifs = [allmotifs[x] for x in idxs]
    else:                             #It's a text string!!
        fname = 'TXT'
        for t in txt.split(','):
            motifs.append(Motif_from_text(t))
    for i in range(len(motifs)): motifs[i].index = i
    for i in range(len(motifs)): motifs[i].file = fname
    return motifs

def pickletxt2motifs(toks):
    """[Utility function] See txt2motifs documentation."""
    fname = toks[0]
    print "# Loading motif pickle from %s"%fname
    F = open(fname,'r')
    DA = pickle.load(F)
    F.close()
    ans = []
    if type(DA) == type({}):
        if len(toks) > 1:
            keys = [x.replace('%',' ') for x in toks[1].split(',')]
            for k in keys: ans.append(DA[k])
        else:
            for k in DA.keys(): DA[k].key = k
            ans = DA.values()
    else: #Assuming DA is a list
        if len(toks) > 1:
            idxs = [int(x) for x in toks[1].split(',')]
            ans  = [DA[x] for x in idxs]
        else:
            ans  = DA
    return ans
    

def sortby(motiflist, property, REV=0):
    """Sort a motif list according to a particular property"""
    mtype = type(Motif())
    for m in motiflist:
        if type(m) != mtype:
            print "Not a Motif Object: ",m
            return
    try:
        motiflist.sort(lambda x,y,p=property: cmp(x.__dict__[p],y.__dict__[p]))
        if REV: motiflist.reverse()
    except:
        print 'Could not sort list.  Probably, the specificied property "%s" is not posessed by all motifs'%property
    

