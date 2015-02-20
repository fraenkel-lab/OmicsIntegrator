===============================
garnet-forest
===============================

Garnet + Forest is a package designed to integrate gene expression
data and/or proteomics data using the protein-protein interaction
network.

Contact: Sara JC Gosline [sgosline@mit.edu], Mandy Kedaigle [mandyjoy@mit.edu]

Copyright (c) 2015 Sara JC Gosline, Mandy Kedaigle
     

System Requirements:
--------------------
1. Python 2.6 or 2.7 (3.x version currently untested): http://www.python.org and the following dependencies (provided with install): 
  - numpy: http://www.numpy.org/
  - scipy: http://www.scipy.org/
  - matplotlib: http://matplotlib.org/
  - Networkx: http://networkx.github.io

2. msgsteiner package: http://areeweb.polito.it/ricerca/cmp/code/bpsteiner

3. Boost C++ library: http://www.boost.org

4. Cytoscape 2.8 or 3.0 for viewing results graphically: http://www.cytoscape.org


Features
--------

* Maps gene expression data to transcription factors using chromatin
  accessibility data

* Identifies proteins in the same pathway as 'hits' using protein
  interaction network

* Integrates numerous high throughput data types to determine testable
  biological hypotheses


Running garnet.py
-----------------
Garnet is a script that runs a series of smaller scripts to map epigenetic data
to genes and then scan the genome to determine the likelihood of a transcription
factor binding the genome near that gene. 

```
Usage: garnet.py [configfilename]

Options:
  -h, --help          show this help message and exit
  --outdir=OUTDIR     Name of directory to place garnet output. DEFAULT:none
  --utilpath=ADDPATH  Destination of chipsequtil library,
                      Default=../src
```

The configuration file should take the following format:

### garnet input

```ini
[chromatinData]
#these files contain epigenetically interesting regions   
bedfile = bedfilecontainingregions.bed   
fastafile = fastafilemappedusinggalaxytools.fasta       
#these two files are provided in the package 
genefile = ../../data/ucsc_hg19_knownGenes.txt  
xreffile = ../../data/ucsc_hg19_kgXref.txt  
#distance to look from transcription start site  
windowsize = 2000  

[motifData]
#motif matrices to be used, data provided with the package
tamo_file = ../../data/matrix_files/vertebrates_clustered_motifs.tamo
#settings for scanning
genome = hg19
numthreads = 4
doNetwork = False
tfDelimiter = .

[expressionData]
expressionFile = tabDelimitedExpressionData.txt
pvalThresh = 0.01
qvalThresh = 
```

#### Chromatin Data

Many BED-formatted (`bedfile`) and FASTA-formatted (`fastafile`) files are included in the examples/ directory. To use your own epigenetic data, upload the BED-file to http://usegalaxy.org and select 
'Fetch Genomic DNA' from the left menu to click on 'Extract Genomic DNA'. This will produce
a FASTA-formatted file that will work with garnet.  We have provided gene (`genefile`) and xref (`xreffile`)  annotations for both hg19 and mm9 - these files can be downloaded from http://genome.ucsc.edu/cgi-bin/hgTables if needed. The `windowsize` parameter determines the maximum distance from a transcription start
site to consider an epigenetic event associated. 2kb is a very conservative metric.

#### motifData

We provide motif data in the proper TAMO format, the user just needs to enter the genome used.
The default `numthreads` is 4, but the user can alter this depending on the processing power 
of their machine. `doNetwork` will create a networkX object mapping transcription factors to 
genes, required input for the [SAMNet algorithm](http://github.com/sgosline/SAMNet).  `tfDelimiter` is an internal parameter to tell garnet how to handle cases when many transcription factors map to the sam
binding motif.

#### expressionData

If the user has expression data to evaluate, provide a tab-delimited file under `expressionFile`. 
File should have two columns, one containing the name of the gene and the second containing the 
log fold change of that gene in a particular condition. We recommend only including those genes 
whose change in expression is statistically significant. P-value (`pvalThresh`) or Q-value 
(`qvalThresh`) thresholds will be used to select only those transcription factors whose 
correlation with expression falls below the provided threshold.

### garnet output
garnet produces a number of intermediate files that enable you to better interpret your data
or re-run a sub-script that may have failed. All files are placed in the directory provided by 
the `--outdir` option of the garnet script.

- **events_to_genes.xls***: This file is a tab-delimited file that takes every epigenetic event and maps it to the closest gene and records its distance to the transcription start site. 

- **events_to_genes_with_motifs.txt**: Result of motif scanning

- **events_to_genes_with_motifs.pkl**: The scanning results are then mapped back to invidual genes
and then compiled into a matrix with each row representing the value for a transcription factor 
binding site and each column representing a value for a gene. This matrix, along with the row 
names (transcritpion factor motif identifies) and column names (genes) are put into a dictionary 
that is then compressed using the Python pickle library.  The individual data are also written out:
  - ***events_to_genes_with_motifs.tgm:*** Numerical values weighting the probability a particular transcription factor is binding to a region near the gene transcription start site.
  - ***events_to_genes_with_motifs_tfids.txt:*** Row names of the .tgm file.
  - ***events_to_genes_with_motifs_geneids.txt:*** Column names of the .tgm file.
- ***events_to_genes_with_motifs_geneid_regression_results.xls***:
- ***events_to_genes_with_motifs_geneid_regression_results_FOREST_INPUT.xls***:



Running forest.py
-----------------
Forest **requires** the msgsteiner package as well as the boost library.

```
Usage: forest.py [options]

Find multiple pathways within an interactome that are altered in a particular
condition using the Prize Collecting Steiner Forest problem

Options:
  -h, --help            show this help message and exit
  -p PRIZEFILE, --prize=PRIZEFILE
                        (Required) Path to the text file containing the
                        prizes. Should be a tab delimited file with lines:
                        "ProteinName PrizeValue"
  -e EDGEFILE, --edge=EDGEFILE
                        (Required) Path to the text file containing the
                        interactome edges. Should be a tab delimited file with
                        3 or 4 columns: "ProteinA        ProteinB
                        Weight(between 0 and 1) Directionality(U or D,
                        optional)"
  -c CONFFILE, --conf=CONFFILE
                        Path to the text file containing the parameters.
                        Should be several lines that looks like:
                        "ParameterName = ParameterValue". Must contain values
                        for w, b, D.  May contain values for optional
                        parameters mu, n, r, g. Default = "./conf.txt"
  -d DUMMYMODE, --dummyMode=DUMMYMODE
                        Tells the program which nodes in the interactome to
                        connect the dummy node to. "terminals"= connect to all
                        terminals, "others"= connect to all nodes except for
                        terminals, "all"= connect to all nodes in the
                        interactome. If you wish you supply your own list of
                        proteins, dummyMode could also be the path to a text
                        file containing a list of proteins (one per line).
                        Default = "terminals"
  --garnet=GARNET       Path to the text file containing the output of the
                        GARNET module regression. Should be a tab delimited
                        file with 2 columns: "TranscriptionFactorName
                        Score". Default = "None"
  --garnetBeta=GB       Parameter for scaling the GARNET module scores. Use to
                        make the GARNET scores on the same scale as the
                        provided scores. Default = 0.01.
  --msgpath=MSGPATH     Full path to the message passing code. Default =
                        "<current directory>/msgsteiner9"
  --outpath=OUTPUTPATH  Path to the directory which will hold the output
                        files. Default = this directory
  --outlabel=OUTPUTLABEL
                        A string to put at the beginning of the names of files
                        output by the program. Default = "result"
  --cyto30              Use this flag if you want the output files to be
                        amenable with Cytoscape v3.0 (this is the default).
  --cyto28              Use this flag if you want the output files to be
                        amenable with Cytoscape v2.8, rather than v3.0.
  --noisyEdges=NOISENUM
                        An integer specifying how many times you would like to
                        add noise to the given edge values and re-run the
                        algorithm. Results of these runs will be merged
                        together and written in files with the word
                        "_noisyEdges_" added to their names. Default = 0
  --shuffledPrizes=SHUFFLENUM
                        An integer specifying how many times you would like to
                        shuffle around the given prizes and re-run the
                        algorithm. Results of these runs will be merged
                        together and written in files with the word

  --knockout=KNOCKOUT   A list specifying protein(s) you would like to "knock
                        out" of the interactome to simulate a knockout
                        experiment, i.e. ['TP53'] or ['TP53', 'EGFR'].
  -k CV, --cv=CV        An integer specifying the k value if you would like to
                        run k-fold cross validation on the prize proteins. 
                        Default = None.
  --cv-reps=CV_REPS     An integer specifying how many runs of cross-
                        validation you would like to run. To use this option,
                        you must also specify a -k or --cv parameter. Default
                        = None.
  -s SEED, --seed=SEED  An integer seed for the pseudo-random number
                        generators. If you want to reproduce exact results,
                        supply the same seed. Default = None.
 
```
                        
### forest input files and parameters

#### Required inputs

The first two options (`-p` and `-e`) are required. You should record your terminal
nodes and prize values in a text file. The file `example/a549/Tgfb_phos.txt` is an example of
what this file should look lie. You should record your interactome and edge
weights in a text file with 3 or 4 columns. We provide the file 
`data/iref_mitab_miscore_2013_08_12_interactome.txt` is an example of this. 

The program will read in these files and create the interactome graph. It will
print warnings whenever it comes across something unexpected, such as an edge
weight that is not between 0 or 1, or a self-edge from one protein to itself.
It will fix these problems and keep going. It will also print a warning if a
large percentage of the names in the prize file do not have matches in the
interactome listed in the edge file. This error may result from using two
different naming schemes for the proteins, or using the wrong interactome for
your purposes.

A sample configuration file, `a549/tgfb_forest.cfg` is supplied. The user can change the
values included in this file directly or can supply their own similarly
formatted file. If the -c option is not included in the command line the
program will attempt to read `conf.txt`. For explanations of the parameters
w (omega), b (beta), and D, see our original publication on this use of the 
PCSF problem. 



#### Optional inputs

The rest of the command line options are optional. 

If you have run the GARNET module to create scores for transcription factors, you can 
include that output file with the `--garnet` option and `--garnetBeta` options. 

The `--dummyMode` option will change which nodes in the terminal are connected to the dummy node in the 
interactome. We provide an example of this using `a549/Tgfb_interactors.txt`. For an explanation of the dummy node, see our original publication on the PCSF problem. 

If the user is not keeping the file `msgsteiner9` in the same directory as forest.py, 
the path needs to be specified using the `--msgpath` option. 

If you would like the output files to be stored in a directory other than the one you are running 
the code from, you can specify this directory with the `--outputpath` option. The names of the 
output files will all start with the word `result` unless you specify another 
word or phrase, such as an identifying label for this experiment or run, with 
the `--outputlabel option`. The `--cyto30` and `--cyto28` tags can be used to 
specify which version of Cytoscape you would like the output files to be 
compatiable with. 

We include two options, `--noisyEdges` and `--shuffledPrizes`, to determine
how robust your results are by comparing them to results with slightly altered 
input values. To use these options, supply a number for either parameter greater than 0. 
If the number you give is more than 1, it will alter values and run the program that number of times and 
merge the results together. The program will add Gaussian noise to the edge 
values you gave in the `-e` option or shuffle the prizes around all the network 
proteins in the `-p` option, according to which option you use. In `--noisyEdges`, the
standard deviation of the Gaussian noise will be the value the user supplied 
for the parameter `n` in the `-c` configuration file, if given. If not given, the standard 
deviation will be the 0.333. The results from these runs will be stored in 
seperate files from the results of the run with the original prize or edge 
values, and both will be outputted by the program to the same directory.

Additional parameters can be set in the configuration file. To reduce the weight 
of hubs we define the parameter `mu`, to assign negative prizes to nodes in the interactome with high degrees (larger mu means a larger penalty for hub nodes), optional parameter `n`, which sets the noise 
level for option `--noisyPrizes`, optional parameter `r`, which
sets the random noise on the edge costs, and optional parameter `g`, which
is a reinforcement parameter that affects convergence.  See the msgsteiner
PNAS publication for details about `r` and `g`.

The knockout option can be used if you would like to simulate a knockout 
experiment by removing a node from your interactome. Specify your knockout 
proteins in a list, i.e. ['TP53'] or ['TP53', 'EGFR'].

The `-k` or `--cv` option can be used if you would like to run k-fold cross 
validation. This will partition the proteins with prizes into k equal 
subsamples. It will run msgsteiner k times, leaving one subsample of prizes out
each time. The `--cv-reps` option can be used if you would like to run k-fold 
cross validation multiple times, each time with a different random 
partitioning of terminals. If you do not supply `--cv-reps` but do provide a k,
cross validation will be run once. Each time it is run, a file called 
`<outputlabel>_cvResults_<rep>.txt` will be created. For each of the k 
iterations, it will display the number of terminals held out of the prizes 
dictionary, the number of those that were recovered in the optimal network as 
Steiner nodes, and the total number of Steiner nodes in the optimal network. 

The `-s` option will supply a seed option to the pseudo-random number 
generators used in noisyPrizes, shuffledPrizes, and the optimization in 
msgsteiner itself. If you want to reproduce exact results, you should supply 
the same seed every time. If you do not supply your own seed, system time is 
used a seed.

###running forest

Once you submit your command to the command line the program will run. It will
display messages as it completes, letting you know where in the process you
are. If there is a warning or an error it will be displayed on the command
line. If the run completes successfully, several files will be created. These
files can be imported into Cytoscape v.3.0 to view the results of the run.
These files will be named first with the outputlabel that you provided (or
"result" by default), and then with a phrase identifying which file type it is.

### forest output

- **objective.txt** contains information about the algorithm run, including any error
messages if there were any during the run.

- **optimalForest.sif** contains the optimal network output of the message-passing
algorithm (without the dummy node). It is a Simple Interaction Format file. To
see the network, open Cytoscape, and click on File > Import > Network >
File..., and then select this file to open. Click OK.

- **augmentedForest.sif** is the same thing, only it includes all the edges in the
interactome that exist between nodes in the optimal Forest, even those edges
not chosen by the algorithm. Betweenness centrality for all nodes was
calculated with this network.

- **dummyForest.sif** is the same as optimalForest.sif, only it includes the dummy
node and all edges connecting to it. This file is useful as a sanity check
(i.e. are there any singleton nodes in your forest, nodes that are only
connected to the network via the dummy node?).

- **edgeattributes.tsv** is a tab-seperated value file containing information for
each edge in the network, such as the weight in the interactome, and the
fraction of optimal networks this edge was contained in (this will be 0 or 1
for a standard run, or something in between if the results are merged together,
i.e. from adding noise to the prizes and re-running the algorithm several
times). To import this information into Cytoscape, first import the network
.sif file you would like to view, and then click on File > Import > Table >
File..., and select this file. Specify that this file contains edge attributes,
rather than node attributes, and that the first row of the file should be
interpreted as column labels. Click OK.

- **nodeattributes.tsv** is a tab-seperated value file containing information for
each node in the network, such as the prize you assigned to it and betweenness
centrality in the augmented network. To import this information into Cytoscape,
first import the network .sif file you would like to view, and then click on
File > Import > Table > File..., and select this file. Specify that this file
contains node attributes, rather than edge attributes, and that the first row
of the file should be interpreted as column labels. Click OK.

When the network and the attributes are imported into Cytoscape, you can alter
the appearance of the network as you usually would using VizMapper.

