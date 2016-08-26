#!/usr/bin/python
'''

GARNET primary script executes 5 sub scripts according to provided configuration file. We provide vertebrate motif data and gene/xref data for hg19 and mm9. 
--------------------------------------------------------------------
Config file:
--------------------------------------------------------------------
Configuration file should provide the following variables.
[chromatinData]
bedfile=[bed file of accessible chromatin regions]
fastafile=[fasta file of same regions, collected via galaxyweb]

genefile=[path to garnet]/examples/ucsc_hg19_knownGenes.txt
xreffile=[path to garnet]/examples/ucsc_hg19_kgXref.txt

windowsize=[distance around transcription starte site]

[motifData]
tamo_file=../data/matrix_files/vertebrates_clustered_motifs.tamo
genome=hg19
numthreads=4
doNetwork=False
tfDelimiter=.

[expressionData]
expressionFile=[name of expression file]
pvalThresh=0.01
qvalThresh=

[regression]
savePlot=False
=======================================================================

'''

__author__='Sara JC Gosline'
__email__='sgosline@mit.edu'

##update this to include direct location of chipsequtil pacakge
import sys,os,re
import argparse
from ConfigParser import ConfigParser

progdir=os.path.dirname(sys.argv[0])


def mapGenesToRegions(genefile,xreffile,bedfile,window='2000',outdir=None):
    '''
    First step of GARNET maps exposed regions in bedfile to closest gene wtihin rpovided window
    calls map_peaks_to_known_genes.py
    '''
    if outdir is None:
        outdir=os.path.splitext(os.path.basename(bedfile))[0]+'eventsWithin'+window
    
    outfile=outdir+'/events_to_genes.xls'
    res=0
    if not os.path.exists(outdir):
        os.system('mkdir '+outdir)

    #old file naming scheme:
    #os.path.splitext(os.path.basename(bedfile))[0]+'eventsWithin'+window+'bp_of_'+os.path.splitext(os.path.basename(genefile))[0]+'.xls'

    ##Step 1: map chromatin regions to nearby genes/transcription start sites
    cmd='python '+os.path.join(progdir,'map_peaks_to_known_genes.py')+' --peaks-format=auto --utilpath='+os.path.join(progdir,'../src/')+' --upstream-window='+window+' --downstream-window='+window+' --tss --map-output='+outfile+' --symbol-xref='+xreffile+' '+genefile+' '+bedfile
    if not os.path.exists(outfile):
        print '\n-----------------------------Gene-region mapping output------------------------------------------\n'
        print 'Running command:\n'+cmd+'\n'

        print 'Mapping genes from '+genefile+' to regions within '+window+' bp of events from '+bedfile+' and putting results in '+outfile
        res=os.system(cmd)
    else:
        print 'File '+outfile+' already exists. If you would like to replace it, delete and re-run'

    return res,outfile


def motifScanning(tamo_file,fastafile,numthreads,genome,closest_gene_file='',gene_list=''):
    '''
    Second step of GARNET scans chromatin regions provided in galaxy-produced FASTA for motif matrix 
    affinity scores
    Arguments:
    tamo_file: TAMO-formatted list of motifs to scan
    fastafile: FASTA-formatted file to scan
    numthreads: number of threads to run, this process can take a while
    genome: which genome build to use
    closest_gene_file: output of map_peaks_to_known_genes.py so that only those fasta sequences that map to genes will be scanned.
    gene_list: list of genes to focus on (i.e. diff ex genes), based on closest_gene_file mappings.
    '''
    if closest_gene_file=='':
        motif_binding_out=re.sub('.fasta','_with_motifs.txt',fastafile)
    else:
        motif_binding_out=re.sub('.xls','_with_motifs.txt',closest_gene_file)


    if os.path.exists(motif_binding_out):
        print '\nIntermediate file '+motif_binding_out+' already exists, if you would like to replace, delete and re-run'
        return 0,motif_binding_out


    scan_cmd='python '+os.path.join(progdir,'motif_fsa_scores.py')+' --motif='+tamo_file+' --genome='+genome+' --outfile='+motif_binding_out+' --genemappingfile='+closest_gene_file+' --scale=10 --threads='+numthreads+' '+fastafile+' --genelist='+gene_list
    print '\n-----------------------------Motif Scanning Output------------------------------------------\n'
    print 'Running command:\n'+scan_cmd+'\n'
    print 'Scanning regions from '+fastafile+' using matrices from '+tamo_file+' and putting results in '+motif_binding_out
    res=os.system(scan_cmd)
    return res,motif_binding_out

def createBindingMatrix(motif_binding_out,outfile,fastafile,tamo_file,use_uniprot=False):
    '''
    Third step of GARNET merges motif scores with closest gene information to create motif/gene 
    scoring matrix with appropriate identifiers
    '''
    if use_uniprot:
        tfs=re.sub('.tamo','_up_tfids.txt',tamo_file)
        matfile=re.sub('.txt','.tgm',motif_binding_out)
    else:
        matfile=re.sub('.txt','.tgm',motif_binding_out)
        tfs=re.sub('.tamo','_tfids.txt',tamo_file)
        
    extra_file=re.sub('.tamo','_source_names.txt',os.path.basename(tamo_file))
    if os.path.exists(extra_file):
        tfs=tfs+','+extra_file
    
    ##using regular gene names here
    map_cmd='python '+os.path.join(progdir,'get_window_binding_matrix.py')+' '+motif_binding_out+' '+outfile+' '+' '+fastafile+" --distance-to-gene='' --motif-id-list="+tfs+' --outfile='+matfile

    pklfile=re.sub('.tgm','.pkl',matfile)
    if os.path.exists(pklfile):
        print '\nIntermediate file '+pklfile+' already exists, if you would like to replace delete and re-run'
        return 0,pklfile

    print '\n-----------------------------Binding Matrix Output------------------------------------------\n'
    print 'Running command:\n'+map_cmd+'\n'
    res=os.system(map_cmd)
    
    return res,pklfile


def getTfsFromRegression(pickle_file,expressionfile,pvalT,qvalT,plot):
    '''
    Fourth step of GARNET is to perform regression with pickled matrix file and expression data
    '''
#    print '\nRunning regression using '+expressionfile+' expression data and '+pickle_file+' binding data'
    outdir=re.sub('.pkl','regression_results.tsv',pickle_file)
#    outdir=os.path.basename(expressionfile).split('.')[-2]+'_'+re.sub('.pkl','',os.path.basename(pickle_file))+'.xls'
    print outdir
    if not os.path.exists(outdir):
        cmd='python '+os.path.join(progdir,'motif_regression.py')+' --outdir='+outdir+' '+pickle_file+' '+expressionfile

        if pvalT is None or pvalT=='':
            if qvalT is None or qvalT=='':
                thresh='0.05'
            else:
                thresh=qvalT
                cmd+=' --use-qval'
        else:
            thresh=pvalT
        cmd+=' --thresh='+thresh
        
        if plot:        
            cmd+=' --plot'

        print '\n-----------------------------Regression Output------------------------------------------\n'
        print 'Running command:\n'+cmd+'\n'
        res=os.system(cmd)
    else:
        res=0
    return res,outdir
    
def main():
    
    srcdir=os.path.join(progdir,'../src')
    
    parser=argparse.ArgumentParser()
    #uniprot option will be deprecated, SAMNet should be able to map to human gene names
    #    parser.add_option('--useUniprot',dest='useUniprot',action='store_true',help='Set this flag to use Uniprot identifies',default=False)
    parser.add_argument('configfilename', help='Path to configuration file.')
    parser.add_argument('--outdir',dest='outdir',help='Name of directory to place garnet output. DEFAULT: none',default=None)
    parser.add_argument('--utilpath',dest='addpath',help='Destination of chipsequtil library, DEFAULT: ../src',default=srcdir)
    parser.add_argument('--allGenes',dest='allgenes',help='Use this flag to use all annotated genes, even if they show no evidence of encoding proteins.',action='store_true',default=False)

    opts=parser.parse_args()
    
    sys.path.insert(0,opts.addpath)
    sys.path.insert(0,opts.addpath+'chipsequtil')

    config=ConfigParser()
    config.read(opts.configfilename)

    ##now check for elements of config file. if they are missing, move onto next step
    ##first step 1 check
    genefile=config.get('chromatinData','genefile')
    bedfile=config.get('chromatinData','bedfile')
    xref=config.get('chromatinData','xreffile')

    window=config.get('chromatinData','windowsize')
    if window is None:
        window='2000'    

    #This variable tracks the results of all the commands. If it becomes non zero, stop.
    keeprunning=0
    
    if genefile is not None and bedfile is not None:
        keeprunning,outfile=mapGenesToRegions(genefile,xref,bedfile,window,opts.outdir)
    else:
        print 'Missing genefile,bedfile or xref file, cannot map genes to regions.'
        sys.exit()
    if keeprunning!=0:
        print 'Error running gene mapping step, check your files and try again'
        sys.exit()
        
    tamofile=config.get('motifData','tamo_file')
    genome=config.get('motifData','genome')
    
    numthreads=config.get('motifData','numthreads')
    if numthreads is None:
        numthreads='1'

    fastafile=config.get('chromatinData','fastafile')

    expr=config.get('expressionData','expressionFile')

    
    ##step 2
    if tamofile is not None and tamofile!='' and genome is not None and fastafile is not None and fastafile!='':
        if os.path.exists(tamofile) and os.path.exists(fastafile):
            keeprunning,binding_out=motifScanning(tamofile,fastafile,numthreads,genome,outfile,expr)
        else:
            binding_out=''
            print 'Missing FASTA file or TAMO file - check your config file and try again.'

        if keeprunning!=0:
            print 'Error running motif-scanning step, check your files and try again'
            sys.exit()
        

    ##step 3
    newfasta=re.sub('.xls','.fsa',outfile)
    if outfile is not None and outfile!='' and binding_out is not None and binding_out!='':
        keeprunning,binding_matrix=createBindingMatrix(binding_out,outfile,newfasta,tamo_file=tamofile,use_uniprot=False)
    else:
        binding_matrix=''

    if keeprunning!=0:
        print 'Error running matrix creation step, check your files and try again'
        sys.exit()

#        pklfile=config.get('motifData','pkl')
    do_network=config.get('motifData','doNetwork')
    delim=config.get('motifData','tfDelimiter')
    if delim is None:##here we want no delimiter if we do not want to tease out individual tfs
        delim=''

    if do_network is not None and do_network!='' and do_network.lower()!='false':
        cmd='python '+os.path.join(progdir,'zipTgms.py')+' --pkl='+binding_matrix+' --genome '+genome+' --as-network --tf-delimiter='+delim
        if opts.allgenes:
            cmd=cmd+' --allGenes'
        print cmd
        os.system(cmd)
        
    pvt=config.get('expressionData','pvalThresh')
    qvt=config.get('expressionData','qvalThresh')
    
    plot = False
    plot_str=config.get('regression','savePlot')
    if plot_str is not None and plot_str != '' and plot_str.lower() != 'false':
        plot = True
    
    ##step 4: regression
    if expr is not None and expr!='':
        #print binding_matrix,expr
        if binding_matrix!='' and os.path.exists(binding_matrix) and os.path.exists(expr):
            keeprunning,tfs=getTfsFromRegression(binding_matrix,expr,pvt,qvt,plot)
        else:
            print 'Cannot perform regression because binding matrix or expression datasets are missing'
    if keeprunning!=0:
        print 'Error running regression step, check your files and try again'
        sys.exit()
   
    
if __name__=='__main__':
    main()
