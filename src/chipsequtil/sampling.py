
import math
import random
import re
import sys
from collections import defaultdict

from chipsequtil import get_org_settings, get_gc_content, get_gc_content_distribution, RefGeneFile
from nib import NibDB, NibException

def kl_divergence(p,q) :
    """Return Kullback-Leibler divergence for two probability distributions
    p and q.  p and q should be indexable objects of the same length where
    p_i corresponds to q_i.
    """
    kl_sum = 0.
    for p_i, q_i in zip(p,q) :
        if p_i != 0 and q_i != 0 :
            kl_sum += p_i * math.log(p_i/q_i)
    return kl_sum

def rejection_sample_bg(fg_dict,organism,bins=100,num_samples=None,verbose=False,
                        bg_match_epsilon=1e-3) :
    '''Generate background sequences according to the size, distance from genes,
    and GC content distributions of the supplied foreground sequences.  *fg_dict*
    is a dictionary of <header>:<sequence> items, where the first part of the
    header must contain:

    >chrX:<start>-<end>

    *organism* is a string that will be used to call the *chipsequtil.get_org
    settings* function and uses the 'genome_dir' and 'annotation_path' keys.
    *bins* is the number of bins to use for representing the GC content
    distribution.  Function returns a dictionary of <header>:<sequence> items
    of generated background sequences.'''

    nib_db = NibDB(nib_dirs=[get_org_settings(organism)['genome_dir']])
    tss_fn = get_org_settings(organism)['annotation_path']
    tss = defaultdict(list)
    for rec in RefGeneFile(tss_fn) :
        tss[rec['chrom']].append((int(rec['txStart']),int(rec['txEnd']),))

    # for each peak find the chromosome, distance to nearest
    # gene, size of peaks in bases, and GC content
    num_samples = len(fg_dict) if not num_samples else num_samples
    dists,sizes=[],[]

    for header,seq in fg_dict.items() :

        # chromosome first field in fasta headers from bed2seq.bedtoseq
        chrom = header.split(':')[0]

        # adjust chromosomes in special cases
        if re.search('random',chrom.lower()) or chrom.lower() == 'chrm' :
            continue

        # start first int in second field of bed2seq.bedtoseq header
        start = int(header.split(':')[1].split('-')[0])
        midpoint = start + len(seq)/2

        # figure out which chromosome we're working on
        tss_chr = tss[chrom]

        # dsts_to_genes is the distance of this peak from all the genes, find minimum
        dists_to_genes = [(s[0]-midpoint) for s in tss_chr]
        try :
            min_dist = min(dists_to_genes,key=lambda x : abs(x))
            dists.append(min_dist)
        except :
            err_str = 'Warning: no genes were found for sequence with header' \
                         ' %s, not using to calculate distributions.\n'%header
            sys.stderr.write(err_str)

        # calculate # bases
        sizes.append(len(seq))

    # GC content distribution for the foreground sequences
    gc_dist = get_gc_content_distribution(fg_dict.values(),bins=bins)

    # max_gc is # peaks w/ highest GC content
    max_gc = max(gc_dist)

    # gene_starts is a list of all genes in (chromosome,gene start) tuples
    gene_starts=[]
    for key in tss.keys():
        chrom=key.split('chr')[-1]
        for x in tss[key]:
            gene_starts.append((key,x[0]))

    # encapsulated function for proposing sequences
    def propose_sequence(dists, gene_starts, sizes, nib_db) :
        # sample a random distance from the list of distances
        d = random.choice(dists)

        # pick a random gene
        chrom, coord = random.choice(gene_starts)

        # propose a starting point for the bg sequence
        midpoint = coord-d+random.randint(-100,100)

        # propose a size for the bg sequence
        size = random.choice(sizes)
        start = int(midpoint-int(size/2))
        stop = int(midpoint+int(size/2))

        #sys.stderr.write("%s:coord=%d size=%d midpoint=%d d=%d\n"%(chrom,coord,size,midpoint,d))
        # if start or stop are negative, skip and try again
        if start < 0 or stop < 0 : seq = None

        # randomly choose strand
        strand = '+' if random.random() > 0.5 else '-'

        # extract the proposed sequence
        try :
            nib_title, seq = nib_db.get_fasta(chrom,start,stop,strand)
        except IOError, e :
            if verbose : sys.stderr.write('IOError in NibDB, skipping: %s,%d-%d,%s\n'%(chrom,start,stop,strand))
            seq = None
        except NibException, e :
            if verbose : sys.stderr.write('NibDB.get_fasta error, %s\n'%e)
            seq = None

        header = '%s:%d-%d'%(chrom,start,stop)

        return header, seq


    # build gc content distribution based on seq length and
    # distance from TSS foreground distributions
    # keep sampling sequences until the distribution stops
    # changing a lot (KL divergence < epsilon)
    bg_gc_cnts = [1.]*bins
    converged = False
    epsilon = bg_match_epsilon
    if verbose : sys.stderr.write('Building empirical background GC content distribution\n')
    while not converged :

        # propose a sequence
        header, seq = propose_sequence(dists,gene_starts,sizes,nib_db)

        # sometimes this happens when there is an error, just try again
        if seq is None :
            continue

        # determine the GC bin for this sequence
        gc_content = get_gc_content(seq)
        gc_bin = -1
        for i in range(bins) :
            win_start = i/float(bins)
            win_end = (i+1)/float(bins)
            if gc_content >= win_start and gc_content < win_end :
                gc_bin = i
                break

        # update the gc content distribution
        sum_cnts = float(sum(bg_gc_cnts))
        if sum_cnts != 0 : # ! on first sequence

            # calculate the current distributions
            last_gc_p = map(lambda x:x/sum_cnts,bg_gc_cnts)
            bg_gc_cnts[gc_bin] += 1
            new_gc_p = map(lambda x:x/sum_cnts,bg_gc_cnts)

            # calculate the kl divergence between last distribution
            # and current one, stopping if less than epsilon
            kl_d = kl_divergence(new_gc_p,last_gc_p)
            if verbose : sys.stderr.write('dist to converge: %.3g\r'%(kl_d-epsilon))
            if kl_d < epsilon :
                converged = True

        else :
            bg_gc_cnts[gc_bin] += 1

    if verbose : sys.stderr.write('\ndone\n')

    # add pseudocounts to account for missing data in bg as to avoid
    # inappropriate scaling in rejection sampling step
    # the fg bin with the largest value that corresponds to an empty
    # bg bin is used to calculate the number of pseudocounts so that
    # the resulting bg bin has the same propotion of counts in it as
    # the original fg bin.  This is calculated as:
    #
    # x_{pseudo} = \frac{p_i\sum_{i=1}^{N}a_i}{1-p_iN}
    #
    # where p_i is the value of the max fg bin w/ zero in the bg bin
    # x_{pseudo} is added to every bin
    pseudocounts = 0
    for fg_i, bg_i in zip(gc_dist,bg_gc_cnts) :
        if fg_i != 0 and bg_i == 0 and fg_i*len(fg_dict) > pseudocounts :
            # if fg_i > 1/sum(bg_gc_cnts) this won't work, but that *shouldn't*
            # ever happen
            if fg_i >= 1./sum(bg_gc_cnts) :
                raise Exception('There was a numeric issue in the rejection sampling routine, please try it again')
            sys.stderr.write(str([fg_i,sum(bg_gc_cnts),len(bg_gc_cnts),1.*fg_i*len(bg_gc_cnts),bg_gc_cnts])+'\n')
            sys.stderr.flush()
            pseudocounts = (fg_i*sum(bg_gc_cnts))/(1-1.*fg_i*len(bg_gc_cnts))

    bg_gc_cnts = map(lambda x: x+pseudocounts/sum(bg_gc_cnts),bg_gc_cnts)
    bg_gc_dist = map(lambda x: x/sum(bg_gc_cnts),bg_gc_cnts)

    # last, find the multiplier that causes the background gc distribution to
    # envelope the foreground gc dist
    z_coeff = gc_dist[0]/bg_gc_dist[0]
    for fg_i, bg_i in zip(gc_dist[1:],bg_gc_dist[1:]) :
        z_coeff = max(z_coeff,fg_i/bg_i)
    bg_gc_dist = map(lambda x: x*z_coeff,bg_gc_dist)

    # start generating bg sequences
    bg_dict = {}

    bg_gcs,bg_sizes=[],[]

    # generate a bg sequence for every fg sequence
    for i in range(num_samples):
        if verbose : sys.stderr.write('%d/%d'%(i,num_samples))

        # propose sequences until one is accepted
        accepted_sequence = False
        while not accepted_sequence:
            if verbose : sys.stderr.write('.')

            # propose a sequence
            header, seq = propose_sequence(dists,gene_starts,sizes,nib_db)

            # problem occured in proposing sequence, just keep going
            if seq is None : continue

            # determine the GC bin for this sequence
            gc_content = get_gc_content(seq)
            gc_bin = -1
            for i in range(bins) :
                win_start = i/float(bins)
                win_end = (i+1)/float(bins)
                if gc_content >= win_start and gc_content < win_end :
                    gc_bin = i
                    continue

            # pick a uniform random number such that it does not exceed
            # the maximum GC content distribution over bins
            # if the random number is <= the GC content for this
            # proposed sequence, accept, otherwise reject
            r = random.random() * bg_gc_dist[gc_bin]
            if r > gc_dist[gc_bin] :
                continue
            else:
                bg_gcs.append(x)
                #bg_sizes.append(size)
                accepted_sequence = True
                bg_dict[header] = seq

        if verbose : sys.stderr.write('\r')
    return bg_dict
