'''
From the original motif scanning, selects those events withing a particular window of a gene and outputs the
tgm, geneids and tf ids for final analysis
'''
__author__='sara jc gosline'
__email__='sgosline@mit.edu'

import re,os,sys
from optparse import OptionParser
from collections import defaultdict
import numpy as np

progdir=os.path.dirname(sys.argv[0])


def build_annotated_tgm(closest_gene_output,distance_to_tss,logistic_score_output,fasta_file,motif_ids,makeWindow=True,tgm_file='',do_pkl=True):
    '''
    Takes existing tgm, and maps to gene names and TF ids within a specific window
    '''
    from chipsequtil import Fasta
    ##get fasta file events, since these are columns in the logistic_score matrix
    seq_ids=Fasta.load(fasta_file,key_func=lambda x: x)

    ##need to get sequence mids in the order they are processed
    ##in the file, this is the index into the score_output file
    ##. ASSUMES GALAXY-formatted FASTA: genome_chr_start_stop_strand[ ignores
    ##things after space]
    seq_mids=[] ##list of FASTA regions, in their appropriate order in the file
    filtered_events={}##gene name of closest gene to event within window
    for k in seq_ids.keys():
        vals=k.split(';')
        if len(vals)==1:
    	    vals=k.split()
        if ':' in vals[0]: #bed tools used
            chr,range=vals[0].split(':')
            low,high=range.split('-')
            mid=str(int(low)+((int(high)-int(low))/2))
            seq_mids.append(chr+':'+mid)
        elif 'random' not in vals[0]: #galaxy tools used
            genome,chr,low,high,strand=vals[0].split('_')
            mid=str(int(low)+((int(high)-int(low))/2))
            seq_mids.append(chr+':'+mid)

        if len(vals)==2:
            filtered_events[chr+':'+mid]=vals[1]
    print 'Found %d events, of which %d have gene names'%(len(seq_mids),len(filtered_events))
    ##this next section relies on xls
    ##filter events that are within distance from closest_gene_output to get gene mapping
    ##
    filtered_fc={}##FC of events within window, in case we want to use in the future

    event_indexes=[] ##

    ##get gene ids, or just use mid of sequence region
    gene_names=[t for t in set(filtered_events.values())]
    print gene_names[0:10]

    #get gene ids for all matrices list loaded in
    mi_files=motif_ids.split(',')
    if len(mi_files)>0:
        #open first motif name file that contains names for each element in TAMO file
        all_tf_names=[a.strip() for a in open(mi_files[0],'rU').readlines()]
    if len(mi_files)>1:
        #if we have additional files, check to see if if names already exist
        for i,f in enumerate(mi_files):
            if i==0:
                next
            try:
                #open file and read in extra ids
                newfs=[a.strip() for a in open(f,'rU').readlines()]
            except:
                print "Error opening file:", sys.exc_info()[0]
                print "Check to make sure file exists at %s"%(f)
                raise

            if len(newfs)==len(all_tf_names):
                #combine existing tf names with these with . delimiter....
                all_tf_names=['.'.join((a,b)) for a,b in zip(all_tf_names,newfs)]

    ##now go through and clean up TF names
    cleaned_tf_names=[]
    for i,a in enumerate(all_tf_names):
        tfn=set([b for b in a.split('.') if '$' not in b and b!=''])
        if(len(tfn)==0):
            tfn=a.split('.')
#        else:
#            print 'Replacing %s with %s'%(a,'.'.join(tfn))
        cleaned_tf_names.append('.'.join(tfn))

    all_tf_names=cleaned_tf_names
    #print len(cleaned_tf_names)


    ##now actually map events to scores
    ##load motif matrix scanning output that maps matrices to regions
    print 'Loading complete motif score file...'
    event_scores=np.loadtxt(logistic_score_output)
    print '\t...Loaded!'

    #create new tgm matrix with approriate file name
    newmat=np.zeros((len(all_tf_names),len(gene_names)),dtype='float')##fill in gene length),dtype='float')
    if makeWindow:
        distance_to_tss=distance_to_tss+'_bpWindow'
    else:
        distance_to_tss=distance_to_tss+'_bpUpstream'

    if tgm_file=='':
        tgm_file=re.sub('.txt','_'+distance_to_tss+'.tgm',os.path.basename(logistic_score_output))
    if do_pkl:
        pkl_file=re.sub('.tgm','.pkl',tgm_file)
    else:
        pkl_file=''

    ##sort event indexes from  that are in the filtered_events file
    event_indexes.sort()

    #populate matrix with greatest score attributed to that gene/tf combo
    for ind,arr in enumerate(event_scores):
        ##name of matrix/motif
        mat=all_tf_names[ind]

        #tfnames=[mat]
        ##here we enumerate which sequences were mapped to a gene within the window
        for k,val in enumerate(seq_mids):#k in event_indexes:

            #here we want the event midpoint for the index
#            val=seq_mids[k]

            #get score for that index
            score=arr[k]

            #now map it to closest gene for that midpoint
            cg=filtered_events[val]

            fc=1.0 ##update this if we want to normalize score by fold change
            score=float(score)*float(fc) ##this should do nothing sine fcgenerally =1

            #if len(tfnames)==1:
            curscore=newmat[all_tf_names.index(mat),gene_names.index(cg)]
            ##updated to include maximum score!!

            if np.abs(score)>np.abs(curscore):
                newmat[all_tf_names.index(mat),gene_names.index(cg)]=score
            #else:
            #    for t in tfnames:
            #        curscore=newmat[all_tf_names.index(t),gene_names.index(cg)]
            #    ##updated to include maximum score!!
            #        if np.abs(float(score))>np.abs(curscore):
            #            newmat[all_tf_names.index(t),gene_names.index(cg)]=float(score)


    ###save these intermediate files for debugging purposes
    np.savetxt(tgm_file,newmat)
    gin=re.sub('.tgm','_geneids.txt',tgm_file)
    tin=re.sub('.tgm','_tfids.txt',tgm_file)

    try:
        open(gin,'w').writelines([g+'\n' for g in gene_names])
        open(tin,'w').writelines([t+'\n' for t in all_tf_names])
    except:
        print "Error opening file:", sys.exc_info()[0]
        print "Check to make sure file exists at %s"%(closest_gene_output)
        raise

    if pkl_file!='':
        zipcmd='python '+os.path.join(progdir,'zipTgms.py')+' '+tgm_file+' '+tin+' '+gin+' --pkl='+pkl_file
        print 'Compressing matrix file into pkl'
        print zipcmd
        os.system(zipcmd)
        return pkl_file
    else:
        return tgm_file

def main():
    '''
    main method
    '''
    usage='usage: %prog [options] motif_scanning_file closest_gene_output fasta_file'
    parser=OptionParser(usage=usage)
    parser.add_option('--distance-to-gene',dest='distance',type='string',default='10000',help='Max distance allowed between event and closest gene')
    parser.add_option('--motif-id-list',dest='motif_ids',type='string',default='',help='Comma-delimited list of files containing TF names that best map to each motif')

    parser.add_option('--utilpath',default=os.path.join(progdir,'../src/'),dest='addpath',help='Destination of chipsequtil library')
    parser.add_option('--outfile',default='',dest='outfile',help='Predefined output file, otherwise will create one automatically')
    parser.add_option('--noPkl',default=True,action='store_false',dest='do_pkl',help='Set this flag if outfile is not in pkl form')

    opts,args=parser.parse_args()

    if len(args)!=3:
        print usage
        exit('Not enough arguments')

    logistic_score_output,closest_gene_output,fasta_file=args

    ##do the path management
    sys.path.insert(0,opts.addpath)
    #sys.path.insert(0,opts.addpath+'chipsequtil')
#    print sys.path
#    print opts.do_pkl

    res=build_annotated_tgm(closest_gene_output,opts.distance,logistic_score_output,fasta_file,opts.motif_ids,tgm_file=opts.outfile,do_pkl=opts.do_pkl)


if __name__=='__main__':
    main()
