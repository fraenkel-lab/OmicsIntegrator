'''
More generic version of the logisitic function

Read in a fasta sequences, scan for every motif in a tamo file, produce numpy matrix as a result

'''

__author__="Sara JC Gosline, Chris W Ng"
__email__="sgosline@mit.edu"

from optparse import OptionParser

#import chipsequtil.Fasta as Fasta
import cPickle
import numpy as np
import math,sys
from multiprocessing import Pool
from csv import DictReader

import os,sys,re

def writeMotNames(m,fname):
    '''
    Write motif names to file
    '''
    ##write out motif names to local file
    newnames=[a.source for a in m]
    genenames=[]
    for n in newnames:
        anns=n.split(';')
        gns=[]
        for i in anns:
            gv=[a.strip() for a in i.split('\t')]
            if len(gv)>1:
                if len(gv)>2 and gv[2] not in gns:
                    gns.append(gv[2])
                elif len(gv)<2:
                    if '$' not in gv[0] and gv[0] not in gns:
                        gns.append(gv[0])
            elif 'Species' not in i and 'included' not in i and i.strip() not in gns:
                gns.append(i.strip())
            #print gns
        gns=set([g.upper() for g in gns if '(' not in g and ')' not in g and ':' not in g and '-' not in g and '/' not in g and 'Delta' not in g and ' ' not in g and '.' not in g])
        genenames.append('.'.join(gns))
    open(fname,'w').writelines([g+'\n' for g in genenames])
    #print fname

#############TRANSFAC SPECIFIC CODE
def load_ids(ids):
    '''
    Loads in TRANSFAC motif identifiers
    '''
    FH=open(ids)
    IDS=[]
    for line in FH:
        mid=line.rstrip('\n')
        IDS.append(mid)
    FH.close()
    return IDS

def motif_matrix(F,motif,outfile,genome,ids,pkl,threads,typ):
    '''
    Creates matrix of motif scan based on TRANSFAC match scores
    '''


    #Load motif and background adjust PSSM
    m=MotifTools.load(motif)
    fname=re.sub('.tamo','_source_names.txt',os.path.basename(motif))
    writeMotNames(m,fname)

#    F=Fasta.load(fsa,key_func=lambda x:x)
    seqs=F.values()
    n_seqs=len(seqs)
    n_motifs=len(m)
    SCORES=np.zeros((n_motifs,n_seqs),dtype='float')
    
    #Load motif ids and profile pickle
    IDS=load_ids(ids)
    PRF=cPickle.load(open(pkl))

    z=zip(m,IDS)
    jobs=[]
    p=Pool(threads)
    for Z in z: jobs.append([Z,PRF,genome,seqs,typ])
    results=p.map(numbs,jobs)
    for i,r in enumerate(results): SCORES[i,:]=r

    np.savetxt(outfile,SCORES,fmt='%.3f')
    
def numbs(args):
    '''
    Calculates the score of all the matches in a particular binding region
    '''
    Z,PRF,genome,seqs,typ=args
    n_seqs=len(seqs)
    M,ID=Z
    thres,SUM,FP=PRF[ID]
    ##if FP is 1, then that can throw things off...
#    if FP==1.0:
#        FP=min(SUM+0.05,0.95)
#        print 'Ajusting FP from 1.0 to '+str(FP)
    ll = M.logP

    if genome in ['hg18','hg19']:
        bg={'A': 0.26005923930888059,
            'C': 0.23994076069111939,
            'G': 0.23994076069111939,
            'T': 0.26005923930888059}
    elif genome in ['mm8','mm9','mm10']: 
        bg={'A': 0.29119881438474354,
        'C': 0.20880118561525646,
        'G': 0.20880118561525646,
        'T': 0.29119881438474354}
    else:
        bg={'A':0.25,'G':0.25,'C':0.25,'T':0.25}

    for pos in ll:
        for letter in pos.keys():
            pos[letter] = pos[letter] - math.log(bg[letter])/math.log(2.0)

    AM = MotifTools.Motif_from_ll(ll)
    mi,ma=AM.minscore,AM.maxscore
    AM.source = M.source
    t=thres*(ma-mi)+mi
    S=np.zeros((1,n_seqs),dtype='float')
    
    #Search every seq for given motif above threshold t and print motif centered results
    for j,seq in enumerate(seqs):
        try:
            seq_fwd = seq.upper()
            matches,endpoints,scores=AM.scan(seq_fwd,threshold=t)
            s=[(x-mi)/(ma-mi) for x in scores]
            aff=affinity(s,SUM,FP,typ)
            #num_bs=len(scores)
            S[0,j]=aff
        except: 
            S[0,j]=0
            #print 'score calc exception',
    return S

def affinity(scores,SUM,FP,typ=6.):
    '''
    takes the afinity based on the FP score (which is just .5?)
    typ is a scaling factor that changes how much to weight the priors
    '''
#    if typ==0:
#        try: w=math.log(9)/(FP-SUM)
#        except: w=math.log(9)/0.1
#        b=math.exp(w*SUM)
#    else: #default to this
    if typ<=1.0:
        typ=1.000001
    try: w=math.log(typ)/(FP-SUM)
    except: w=math.log(typ)/0.1
    #SJCG: altered prob of being unbound to scale with the typ
    ##scaling with length of scores penalizes lower probability sites
    ##    b=math.log(10.0-typ)*math.exp(w*SUM)*len(scores)
    #b=2*(10-typ)*math.exp(w*SUM) this looks pretty good, lowers total scores
    b=2*typ*math.exp(w*SUM)
        #try: w=2*math.log(7/3.)/(FP-SUM)
        #except: w=2*math.log(7/3.)/0.1
        #b=7/3*math.exp(w*SUM)
    a=0
    for x in scores: a+=math.exp(w*x) #probability of being bound at any region
    A=a/(b+a)#prob of being bound (a) / prob unbond (b) + prob bound (a)
    return A

def reduce_fasta(fsa_dict,gene_file,gene_list):
    '''
    Takes FASTA file and reduces events to those found in gene_file
    '''
    #read in genes from differential expression data
    include_genes=[]
    if os.path.exists(gene_list):
        include_genes=[a.strip().split()[0] for a in open(gene_list,'rU').readlines()]
    else:
        sys.exit('ERROR: Cannot find file %s'%gene_list)
        
    if len(include_genes)>0:
        print 'Found %d genes in differential expression data, reducing FASTA to only include regions nearby'%(len(include_genes))
    else:
        sys.exit('ERROR: Found 0 differentially expressed genes in %s, check file format.'%gene_list)

    #read in genes mapped to peaks
    if os.path.exists(gene_file):
        closest_gene=DictReader(open(gene_file,'rU'),delimiter='\t')
    else:
        sys.exit('ERROR: Cannot find expected output file from step 1, %s'%gene_file)
        
    #get midpoints of genes near peaks and differentially expressed
    mapped_mids=set()
    mid_to_gene={}
    decount=0
    totcount=0
    for g in closest_gene:
        totcount += 1
        if g['geneSymbol'] not in include_genes:
           # gene near peak, but no expression data
           continue
        try: #this will work for bed
            midval = g['chrom']+':'+str(int(g['chromStart'])+(int(g['chromEnd'])-int(g['chromStart']))/2)
        except:
            try: #this will work for MACS
                midval = g['chr']+':'+str(int(g['start'])+(int(g['end'])-int(g['start']))/2)
            except:#this will work for GPS
                midval = 'chr'+g['Position']
        mapped_mids.add(midval)
        if 'geneSymbol' in g.keys():
            mid_to_gene[midval] = g['geneSymbol']
        else:
            mid_to_gene[midval] = g['knownGeneID']
        decount=decount+1

    if totcount == 0:
        sys.exit('ERROR: Output file from step 1, %s, is empty or formatted incorrectly.'%gene_file)
    if decount==0:
        sys.exit('There were zero differentially expressed genes found near the epigenetic regions in the windowsize you provided.')

    #reduce fasta dict to only get those genes that are differentially expressed
    new_seq={}    
    for k in fsa_dict.keys():
        vals=k.split(';')
        if len(vals)==1:
            vals=k.split(' ')
        #find midpoint of regions in fasta file
        if ':' in vals[0]: #just in case bedtools were used 
            chr,range=vals[0].split(':')
            low,high=range.split('-')
            mid=str(int(low)+((int(high)-int(low))/2))
            seq_mid=chr+':'+mid
        elif 'random' not in vals[0]: #galaxy tools used, ignoring random chromosomes
            allvals=vals[0].split('_')
            if(len(allvals)<4):
                print 'Cannot find sequence data for '+vals[0]
            else:
                genome=allvals[0]
                chr=allvals[1]
                low=allvals[2]
                high=allvals[3]
            mid=str(int(low)+((int(high)-int(low))/2))
            seq_mid=chr+':'+mid
        if seq_mid in mapped_mids:
            #new_seq will hold regions and gene names mapped to sequence
            new_seq[k+' '+mid_to_gene[seq_mid]]=fsa_dict[k] 

    print 'Found '+str(len(new_seq))+' events from FASTA file that map to '+str(decount)+' event-gene matches  and '+str(len(include_genes))+' genes out of '+str(len(fsa_dict))+' events'

    if len(new_seq) == 0:
        sys.exit('ERROR: There were zero sequences found in the FASTA file that map to the right regions. The FASTA file should contain the sequences for the regions in your epigenetic file.')

    ##now write new fasta file
    Fasta.write(new_seq,re.sub('.xls','.fsa',gene_file))
    return new_seq
    

##get program directory
progdir=os.path.dirname(sys.argv[0])

def main():
    usage = "usage: %prog [opts] fasta_file"
    srcdir=os.path.join(progdir,'../src')

    parser=OptionParser(usage)
    
    parser.add_option("--motif", dest="motif",default=os.path.join(progdir,"../data/matrix_files/vertebrates_clustered_motifs.tamo"),help='The .tamo formatted motif file to use for motif matching and scoring')
    parser.add_option('--scores',dest='pkl',default=os.path.join(progdir,'../data/matrix_files/motif_thresholds.pkl'),help='PKL file of matrix score thresholds')
    parser.add_option('--ids',dest='ids',default=os.path.join(progdir,'../data/matrix_files/vertebrates_clustered_motifs_mIDs.txt'),help='List of Exemplar motifs in motif cluster')
    
    parser.add_option('--genemappingfile',dest='gene_file',default='',help='File indicating which regions are mapped to genes, enabling the reduction of the FASTA file for gene-relevant regions')
    parser.add_option("--genome", dest="genome", default='mm9',help='The genome build that you are using, used to estimate binding site priors')
    parser.add_option('--utilpath',dest='addpath',default=srcdir,help='Destination of chipsequtil library, Default=%default')
    parser.add_option('--genelist',dest='genelist',default='',help='List of genes (will select first column if multiple) to include in scan based on --genemappingfile')
    parser.add_option("--outfile", dest="outfile")

#    parser.add_option('--logistic',dest='logistic',action='store_true',default=False,help='Set to true to scale multiple matches into a logistic curve')
    parser.add_option('--threads',dest='threads',type='string',default='4',help='Set number of threads if using logistic scoring')
    parser.add_option('--scale',dest='typ',type='string',default='6')
    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
        
    fsa=args[0]
    motiffile=opts.motif
    ##append path to chipsequtil/TAMO
    sys.path.insert(0,opts.addpath)
    global MotifTools
    from chipsequtil import motiftools as MotifTools
    global Fasta
    from chipsequtil import Fasta
#    sys.path.insert(0,opts.addpath+'chipsequtil')
    
    fsa_dict=Fasta.load(fsa,key_func=lambda x:x)
    if opts.gene_file!='':
        print 'Reducing FASTA file to only contain sequences from '+opts.gene_file
        fsa_dict=reduce_fasta(fsa_dict,opts.gene_file,opts.genelist)


    motif_matrix(fsa_dict,motiffile,opts.outfile,opts.genome,opts.ids,opts.pkl,int(opts.threads),typ=float(opts.typ))

if __name__=='__main__':
    main()
