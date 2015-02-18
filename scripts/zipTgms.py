'''
Quick utility script that takes three files and zips them into a single pkl file
'''

import sys,re,pickle
import numpy,os
from optparse import OptionParser
import networkx
from csv import DictReader

progdir=os.path.dirname(sys.argv[0])
def get_transcriptional_network_from_tgm(tgm,addmrna=True,score_thresh=0.1,expressed_genes=set(),tf_annotation=set()):

    # load in values from TGM
    scores=tgm['matrix']
    tfs=tgm['tfs']
    gs=tgm['genes']
    delim=tgm['delim']

    # combine tf and their scores
    tf_scores = zip(tfs,scores)



    ##now initialize network
    transcription_graph = networkx.DiGraph()
    all_edge_weights=[] ##keep these for normalization
    ##remove this
    p300counts=0
    p300ids=['EP300','Ep300','icrogid:820620', 'icrogid:4116836','EP300_HUMAN']
    print "Creating network digraph for "+str(len(tfs))+' matrices and '+str(len(gs))+' genes from '+str(len(expressed_genes))+' expressed genes'
    count=0
    # Loop through and first expand tf set based on multiple tfs in the same family, ex: EGR1.EGR2
    for tfs,score in tf_scores:
        count+=1
        if count%50==0:
            print 'Processed '+str(count)+' of '+str(len(tf_scores))+' TF motif clusters.'
            print '...tf network has '+str(transcription_graph.number_of_edges())+' edges'
        
        if delim!='' and delim in tfs: # In the case of a multi-tf family, break them up
#            print 'Splitting tfs by '+delim
            alltfs=tfs.split(delim)
        else:
#            print 'Not splitting tfs: '+tfs
            alltfs=[tfs]
        for tf in alltfs:
            if tf in p300ids or tf.upper() in p300ids:
                p300counts+=1
                continue
            if len(tf_annotation)>0 and tf not in tf_annotation:
                ##removed for now, dont want to filter TFs, want to include all
                continue
            gscore=zip(gs,score)##zip each gene name to score vector
            for g,sc in gscore:
                if len(expressed_genes)>0 and g not in expressed_genes:
                    continue
                if sc<=score_thresh:
                    continue
                if addmrna:
                    g+='mrna'
                if sc==1.0:
                    sc=0.9999
                all_edge_weights.append(sc)
                #now we can finally add the interaction
                try:
                    w=transcription_graph[tf][g]['weight']
                except KeyError:
                    transcription_graph.add_edge(tf,g,{'weight':sc})
                    w=sc
                if sc > w:
                    transcription_graph[tf][g]['weight']=sc


    print 'Removed '+str(p300counts)+' interactions containing p300'
        
    print 'Returning transcriptional network with '+str(transcription_graph.number_of_nodes())+' nodes and '+str(transcription_graph.number_of_edges())+' edges'
    return transcription_graph


if __name__=='__main__':
    usage='USAGE: python zipTgms.py [tgmfile] [tfids] [geneids]\nFile arguments can be left blank if --pkl option contains pickled dictionary containing files (from get_window_binding_matrix.py)'
    parser=OptionParser(usage=usage)
    parser.add_option('--pkl',dest='pkl',type='string',help='Name of pkl file to store combined object in, or name of existing pkl if using just network option.')
    parser.add_option('--tf-delimiter',dest='delim',type='string',help='Delimiter used to separate TF names. If left blank, will treat matrix as TF and create network with matris as node. DEFAULT: \'.\'')
    parser.add_option('--as-network',dest='as_network',action='store_true',default=False,help='Set this flag to save network as networkx object instead of matrix/names')
    parser.add_option('--genome',dest='genome',default='hg19',help='Genome build for species-specific filtering')
    parser.add_option('--minscore',dest='minscore',default='0.3',type='string',help='If building networkX object, will remove edges with weight less than this value')
    opts,args=parser.parse_args()

    if len(args)!=3:
        if os.path.exists(opts.pkl) and opts.as_network:
            print 'Found existing pkl '+opts.pkl+', making into networkX object'
            resfile=pickle.load(open(opts.pkl,'rU'))
            mat=resfile['matrix']
            tfs=resfile['tfs']
            geneids=resfile['genes']
            tf_delimiter=resfile['delim']
            if opts.delim is not None:
                tf_delimiter=opts.delim
                print 'Over-riding existing delimiter ('+tf_delimiter+') with new one: '+opts.delim
            fname=re.sub('pkl','thresh'+opts.minscore+'network.pkl',opts.pkl)
        else:
            print usage
            sys.exit('Need 3 arguments or pkl containing files in dictionary')
    else:
        tgm,tfids,geneids=args
        print 'Loading files...'
        print 'TGM: '+tgm
        ##first load tgm into numpy matrix
        mat=numpy.loadtxt(tgm)
        print 'TFs: '+tfids
        tfs=[a.strip() for a in open(tfids,'rU').readlines()]
        print 'Genes: '+geneids
        geneids=[a.strip() for a in open(geneids,'rU').readlines()]
        if opts.delim is None:
            tf_delimiter='.'
        else:
            tf_delimiter=opts.delim

        fname=opts.pkl


    ##create dictionary object for regression
    resfile={'matrix':mat,'tfs':tfs,'genes':geneids,'delim':tf_delimiter}


    if opts.as_network: ##now create network DiGraph if desired
        #first filter by expressed proteins
        prots={'mouse':set(),'human':set(),'mammal':set()}
        res=DictReader(open(progdir+'/../data/mouse-human-uniprot2014-1.tab','rU'),delimiter='\t')
        ##just add in all identifiers
        for row in res:
            spec=row['Organism']
            ps=set([row['Entry'],row['Entry name']])|set(row['Gene names'].split())
            prots['mammal']|=set(ps)
            if spec=='Mus musculus (Mouse)':
                prots['mouse']|=set(ps)
            else:
                prots['human']|=set(ps)

        if opts.genome in ['hg19','hg18','hg17']:
            spec='human'
        elif opts.genome in ['mm8','mm9','mm10']:
            spec='mouse'
        
            ##now get species list of proteins
        spec_prot=[]
        if tf_delimiter!='':
            spec_prot=prots['human']

        resfile=get_transcriptional_network_from_tgm(resfile,score_thresh=float(opts.minscore),expressed_genes=prots[spec],tf_annotation=spec_prot)

    pickle.dump(resfile,open(fname,'w'))
    print 'Combined file saved to '+fname
