'''
File to handle motif/expression regression
'''

__author__='Anthony Soltis'
__email__='asoltis@mit.edu'

import sys,os,pickle,re
from optparse import OptionParser
import numpy as np
from scipy import stats
import fileinput
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt

def load_tgm(tgm_fn):
    '''
    Load tgm file and produce output matrix.
    Output is transposed numpy array object.
    '''
    print 'Loading tgm file...'
    tgm = []
    for line in fileinput.input(tgm_fn):
        l = line.strip('\n').split()
        tgm.append(l)

    # display results, return array
    s = np.asarray(tgm).T.shape
    print 'TGM file loaded with %d genes by %d motifs.'%(s[0],s[1])
    return np.asarray(tgm).T

def load_ids(ids_fn):
    '''
    Load ids filename and store as list.
    '''
    ids = []
    for line in fileinput.input(ids_fn):
        l = line.strip('\n')
        ids.append(l)

    return ids

def load_response(data_fn):
    '''
    Load ydata and return numpy vector.
    Input file should have one value per-row.
    '''
    r_data = []
    r_genes = []
    for line in fileinput.input(data_fn):
        row=line.strip('\n').split('\t')
        if len(row)>1:
            r_genes.append(row[0])
            r_data.append(float(row[1]))
        else:
            r_data.append(float(row[0]))
#        r_data.append(float(line.strip('\n')))
    
    print 'Response data file loaded with %d values.'%(len(r_data))
    return np.asarray(r_data),r_genes

def map_data(Xdata,Xnames,Ydata,Ynames):
    '''
    Map X (predictor) data to Y (response) data using X and Y data ids (i.e. gene names).
    '''
    # Intersect two gene lists
    Xinds = []
    Yinds = []
    for i,Xgene in enumerate(Xnames):
        for j,Ygene in enumerate(Ynames):
            if Xgene == Ygene:
                Xinds.append(i)
                Yinds.append(j)
    Xdata_out = Xdata[Xinds,:]
    Ydata_out = Ydata[Yinds]
    return Xdata_out,Ydata_out

def perform_regression(X,Y,motif_ids,norm):
    '''

    '''
    reg_results = []
    for i in range(0,X.shape[1]):
        # Set up data
        x = np.array(X[:,i],dtype=float)
        if norm != None:
            if norm == 'log2':
                y = np.log2(Y+.1)
            elif norm == 'log10':
                y = np.log10(Y+.1)
        else: y = Y

        # Perform regression
        slope,intercept,r_val,p_val,std_err = stats.linregress(x,y)
        reg_results.append(([motif_ids[i],slope,p_val,i]))

        #fig = plt.figure()
        #ax1 = fig.add_subplot(111)
        #ax1.plot(x,y,'bo',x,intercept+slope*x,'k')
        #ax1.set_title(motif_ids[i])
        #fig.savefig(open('test_plots/'+motif_ids[i]+'.pdf','w'),dpi=300)
        #plt.close()

    return sorted(reg_results,key=lambda x: x[2])

def fdr_correction(results):
    '''
    Compute FDR corrected p-values based on Benjamini-Hochberg procedure.
    '''
    new_results = []
    num_tests = len([r for r in results if str(r[1])!='nan'])
    print 'Correcting for '+str(num_tests)+' numeric values'
    for i in range(0,num_tests):
        tup = results[i]
        pval = tup[2]
        fdr = num_tests*pval/(i+1)
        if fdr > 1.0: fdr = 1.0
        tup+=(fdr,)
        new_results.append(tup)

    return new_results

def main():
    
    usage = "%prog [options] <scores.tgm or scores.tgm.pkl> <response_values.tab>"
    description = "Script that takes a predicted TF-Gene matrix and uses a linear regression to identify which TFs have binding scores correlated with gene expression changes."
    parser = OptionParser(usage=usage,description=description)

    ##get program directory
    progdir=os.path.dirname(os.path.abspath(sys.argv[0]))
    
    # Options
    parser.add_option('--outdir','--out',dest="outdir",default='./test_out.txt',
                      help='Choose output file name. Default is %default.')
    parser.add_option('--motif-ids','--motif-ids',dest='motif_ids',default=None,
                      help='OPTIONAL: If input file is in text format (.tgm), provide motif ids corresponding to tgm file motifs.')
    parser.add_option('--tgm-genes',dest='tgm_genes',default=None,
                      help='OPTIONAL: If input file is in text format (.tgm), provide gene ids corresponding to tgm file genes.')
    parser.add_option('--response-genes',dest='response_genes',default=None,
                      help='OPTIONAL: If two-column file is not provided, add in gene ids corresponding to response values.')
    parser.add_option('--norm-type',dest='norm_type',default=None,
                      help='Choose normalization type for response data. Choices are: "log2", "log10".\
                            Default is %default.')    
    parser.add_option('--use-qval',dest='use_qval',action='store_true',default=False,help='If set this the FOREST input file will contain -log(qval) instead of -log(pval). Default:%default')
    parser.add_option('--thresh',dest='thresh',type='string',default='0.9',help='P/Q-Value threshold to illustrate results. Default:%default')
    parser.add_option('--gifdir',dest='motifs',default=os.path.join(progdir,'../data/matrix_files/gifs'),
                      help='Directory containing motif GIFs to illustrate results. Default is %default')

    # get options, arguments
    (opts,args) = parser.parse_args()

    # Handle arguments
    tgm_fn = args[0]
    response_data_fn = args[1]

    # Load in Y-vector data (gene expression, fold-changes, etc.)
    response_data,response_genes = load_response(response_data_fn)

    print 'Trying to get file type...'
    ext=tgm_fn.split('.')[-1]
    if ext.lower()=='pkl':
        print '...found PKL file'
        pkl=True
    else:
        print '...found text file, looking for additional data files in options'
        pkl=False

    # Handle options
    outdir = opts.outdir
    motif_ids = opts.motif_ids
    if motif_ids == None and not pkl:
        print 'Must provide motif ids file or use pickled dictionary. Exiting.'
        sys.exit()
    tgm_genes = opts.tgm_genes
    if tgm_genes == None and not pkl:
        print 'Must provide gene ids for motifs file or use pickled dictionary. Exiting.'
        sys.exit()
#    response_genes = opts.response_genes
    if opts.response_genes == None and len(response_genes)==0:
        print 'Must provide gene ids for response data or have a two-column data file. Exiting.'
        sys.exit()
    norm_type = opts.norm_type
    valid_norm_types = ['log2','log10']
    if norm_type != None:
        if norm_type not in valid_norm_types:
            print 'Normalization type not valid. Exiting.'
            sys.exit()
    
    if pkl:
        #load in values from dictionary
        tgmdict=pickle.load(open(tgm_fn,'rU'))
        tgm_data=tgmdict['matrix'].T
        motif_ids=tgmdict['tfs']
        tgm_genes=tgmdict['genes']
        delim=tgmdict['delim']
    else:
        # Load in transcription factor affinity matrix and IDs
        tgm_data = load_tgm(tgm_fn)
        motif_ids = load_ids(motif_ids)
        tgm_genes = load_ids(tgm_genes)
        delim='.'

    #now load response_genes if they're not loaded yet
    if len(response_genes)==0:
        response_genes = load_ids(opts.response_genes)


    # Map predictor data to response data
    X,Y=map_data(tgm_data,tgm_genes,response_data,response_genes)
    
    # Perform regression
    reg_results=perform_regression(X,Y,motif_ids,norm_type)
    
    # FDR correction
    new_results = fdr_correction(reg_results)


    dn=os.path.dirname(outdir)
    if dn!='' and dn!='./' and not os.path.exists(dn):
        os.system('mkdir '+dn)
        
    # Write to TEXT file complete results
    of = open(outdir,'w')
    of.writelines('\t'.join(['Motif','Slope','p-val','q-val'])+'\n')
    for res in new_results:
        if str(res[1])=='nan':
            continue    
        ostr = '\t'.join([res[0],str(res[1]),str(res[2]),str(res[4])]) + '\n'
        of.writelines(ostr)
    of.close()

    ##now create HTML writeup
    of= open(re.sub(outdir.split('.')[-1],'html',outdir),'w')
    of.writelines("""<html>
                     <title>GARNET Results</title>
                     <h3>GARNET regression results</h3> 
                     <p>This table includes the results for GARNET TF-motif discovery and regression. This Table includes the non-zero results of the linear regression</p>
                     <table width="90%">
                      <tr><th style="width:25%">Motif Cluster</th><th style="width:12%">Slope</th><th style="width:12%">P-value</th><th style="width:12%">Q-value</th><th style="width:35%">LOGO</th></tr>
                """)
    for res in new_results:
        if str(res[1])=='nan':
            continue    
        if res[4]<float(opts.thresh):
            motifgif=os.path.join(opts.motifs,'motif'+str(res[3])+'.gif')
            ostr = "<tr><td>"+' '.join(res[0].split('.'))+"</td><td>"+str(res[1])+'</td><td>'+str(res[2])+"</td><td>"+str(res[4])+"</td><td><img src=\""+motifgif+"\" scale=80%></td></tr>\n"
            of.writelines(ostr)
    of.writelines("</table></html>")
    of.close()
    
    
    ##now write to FOREST-friendly input file
    ##collect dictionary of all individual tf names and their regression p-values
    regdict={}
    for row in new_results:
        tfs=[t for t in row[0].split(delim) if t!='' and ' ' not in t]
	#print row
        if str(row[1])=='nan':
            continue    
        if opts.use_qval:
            if row[4]>float(opts.thresh):
                continue
        elif row[2]>float(opts.thresh):
            continue
        for tf in tfs:
            if row[2]==1:
                continue
            if opts.use_qval:
                lpv=-1.0*np.log2(float(row[4]))#calculate neg log pvalue
            else:
                lpv=-1.0*np.log2(float(row[2]))#calculate neg log pvalue
            try:
                cpv=regdict[tf]
            except KeyError:
                cpv=0.0
            if lpv>cpv:
                regdict[tf]=lpv
    print 'Found '+str(len(regdict))+'Tf scores for '+str(len(new_results))+' motif results'
    of=open('FOREST_input_'+outdir,'w')
    for tf in regdict.keys():
        val=regdict[tf]
        of.write(tf+'\t'+str(val)+'\n')
    of.close()
        
if __name__ == '__main__': main()
