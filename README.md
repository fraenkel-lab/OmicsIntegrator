<center><img src="http://fraenkel-nsf.csbi.mit.edu/omicsintegrator/omicsI_logo.png" height="40%" width="40%" ></center>

===============================
Omics Integrator
===============================

[![Build Status](https://travis-ci.org/fraenkel-lab/OmicsIntegrator.svg?branch=master)](https://travis-ci.org/fraenkel-lab/OmicsIntegrator)

Omics Integrator is a package designed to integrate proteomic data, gene expression data and/or epigenetic data using a protein-protein interaction network. It is comprised of two modules, Garnet and Forest.

Contact: Amanda Kedaigle [mandyjoy@mit.edu]

Copyright (c) 2015 Massachusetts Institute of Technology
All rights reserved.

Reference:
--------------------
[Network-Based Interpretation of Diverse High-Throughput Datasets through the Omics Integrator Software Package](http://dx.doi.org/10.1371/journal.pcbi.1004879)
Tuncbag N<sup>\*</sup>, Gosline SJC<sup>\*</sup>, Kedaigle A, Soltis AR, Gitter A, Fraenkel E. *PLoS Comput Biol* 12(4): e1004879. doi:10.1371/journal.pcbi.1004879.

System Requirements:
--------------------
1. Python 2.6 or 2.7 (3.x version currently untested) and the dependencies
below. We recommend that users without an existing Python environment
install Anaconda (https://www.continuum.io/downloads) to obtain Python
2.7 and the following required packages:
  - numpy: http://www.numpy.org/
  - scipy: http://www.scipy.org/
  - matplotlib: http://matplotlib.org/
  - Networkx: http://networkx.github.io

2. msgsteiner package (version 1.3): [code](http://staff.polito.it/alfredo.braunstein/code/msgsteiner-1.3.tgz), [license](http://areeweb.polito.it/ricerca/cmp/code/bpsteiner)

3. Boost C++ library: http://www.boost.org

4. Cytoscape for viewing results graphically (tested on versions 2.8-3.2):
http://www.cytoscape.org

Features
--------

* Maps gene expression data to transcription factors using chromatin
  accessibility data

* Identifies proteins in the same pathway as `hits` using protein interaction
  network

* Integrates numerous high throughput data types to determine testable
  biological hypotheses

Installation:
--------------------
Omics Integrator is a collection of Python scripts and data files so can be
easily installed on any system. Steps 1 through 4 are only required for Forest,
and you may skip to step 5 if you will only be running Garnet.

1. Boost is pre-installed on many Linux distributions. If your operating system
does not include Boost, follow the [Boost getting started
guide](http://www.boost.org/doc/libs/1_59_0/more/getting_started/index.html) for
instructions on how to download the library and extract files from the archive.
To use the [Homebrew](http://brew.sh/) package manager for Mac simply type `brew install boost` to install the library.
2. Download `msgsteiner-1.3.tgz` from http://staff.polito.it/alfredo.braunstein/code/msgsteiner-1.3.tgz ([license](http://areeweb.polito.it/ricerca/cmp/code/bpsteiner))
3. Unpack files from the archive: `tar -xvf msgsteiner-1.3.tgz`
4. Enter the `msgsteiner-1.3` subdirectory and run `make`
  * See [this advice](./patches) on compiling the C++ code if you encounter problems and [this advice](https://github.com/fraenkel-lab/OmicsIntegrator/issues/22) regarding compilation issues on OS X.
  * Make a note of the path to the compiled msgsteiner file that was created, which you will use when running Forest.
  * In Linux, use `readlink -f msgsteiner` in the `msgsteiner-1.3` subdirectory to obtain the path.
5. Download the Omics Integrator package: [OmicsIntegrator-0.3.0.tar.gz](./dist/OmicsIntegrator-0.3.0.tar.gz)
6. Unpack files from the archive: `tar -xvzf OmicsIntegrator-0.3.0.tar.gz`
7. Make sure you have all the requirements using the pip tool by entering the
directory and typing: `pip install -r requirements.txt`
  * Some users have reported errors when using this command to install matplotlib. To fix, install matplotlib independently (http://matplotlib.org) or use Anaconda as indicated above.

Now Omics Integrator is installed on your computer and can be used to analyze
your data.

Examples
-----------------
We provide many scripts and files to showcase the various capabilities of Omics
Integrator.  To run this:

1. Download the [example files](./dist/OmicsIntegratorExamples.tar.gz)
2. Unpack by typing `tar -xvzf OmicsIntegratorExamples.tar.gz` in the `dist`
directory.
3. Move the unpacked files into the `example` directory.

For specific details about the examples, check out the [README
file](./example/README.md) in the example directory.

Running garnet.py
-----------------

Garnet is a script that runs a series of
smaller scripts to map epigenetic data to genes and then scan the genome to
determine the likelihood of a transcription factor binding the genome near that
gene.

```
Usage: garnet.py [configfilename]

  -s SEED, --seed=SEED  An integer seed for the pseudo-random number
                        generators. If you want to reproduce exact results,
                        supply the same seed. Default = None.


Options:
  -h, --help            show this help message and exit
  --outdir=OUTDIR       Name of directory to place garnet output. DEFAULT:none
  --utilpath=ADDPATH    Destination of chipsequtil library, Default=../src
```

Unlike Forest, the Garnet configuration file is a positional argument and must not
be preceded with `--conf=`.  The configuration file should take the following format:

### garnet input

```
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

[regression]
#for generating and saving regression plots
savePlot=False
```

#### Chromatin Data

Many BED-formatted (`bedfile`) and FASTA-formatted (`fastafile`) files are
included in the examples/ directory. `bedfile` can also be output from MACS
(with a `.xls` extension) or GPS/GEM (with a `.txt` extension).
To use your own epigenetic data, convert to BED and upload the
BED-file to http://usegalaxy.org and select `Fetch Alignments/Sequences` from the left
menu to click on `Extract Genomic DNA`. This will produce a FASTA-formatted file
that will work with garnet.  We have provided gene (`genefile`) and xref
(`xreffile`)  annotations for both hg19 and mm9 - these files can be downloaded
from http://genome.ucsc.edu/cgi-bin/hgTables if needed. The `windowsize`
parameter determines the maximum distance from a transcription start site to
consider an epigenetic event associated. 2kb is a very conservative metric.

#### motifData

We provide motif data in the proper TAMO format, the user just needs to enter
the genome used.  The default `numthreads` is 4, but the user can alter this
depending on the processing power of their machine. `doNetwork` will create a
NetworkX object mapping transcription factors to genes, required input for the
[SAMNet algorithm](http://github.com/sgosline/SAMNet).  `tfDelimiter` is an
internal parameter to tell Garnet how to handle cases when many transcription
factors map to the sam binding motif.

#### expressionData

If the user has expression data to evaluate, provide a tab-delimited file under
`expressionFile`.  File should have two columns, one containing the name of the
gene and the second containing the log fold change of that gene in a particular
condition. We recommend only including those genes whose change in expression is
statistically significant. P-value (`pvalThresh`) or Q-value (`qvalThresh`)
thresholds will be used to select only those transcription factors whose
correlation with expression falls below the provided threshold.

#### regression

Linear regression plots are placed in a subdirectory named `regression_plots` if
`savePlot=True` in the configuration file.

### Garnet output

Garnet produces a number of intermediate files that enable you
to better interpret your data or re-run a sub-script that may have failed. All
files are placed in the directory provided by the `--outdir` option of the
garnet script.

- **events_to_genes.fsa**: This file contains the regions of the fastafile
  provided in the configuration file that are within the specified distance to a
  transcription start site.

- **events_to_genes.xls**: This file contains each region, the epigenetic
  activity in that region, and the relationship of that region to the closest
  gene.

- **events_to_genes_with_motifs.txt**: This contains the raw transcription
  factor scoring data for each region in the fasta file.

- **events_to_genes_with_motifs.tgm**: This contains the transcription factor
  binding matrix scoring data mapped to the closest gene.

- **events_To_genes_with_motifs_tfids.txt**: Names of transcription factors (or
  columns) of the matrix.

- **events_to_genes_with_motifs_geneids.txt**: Names of genes (or rows) of the
  matrix.

- **events_to_genes_with_motifs.pkl**: A Pickle-compressed Python File
  containing a dictionary data structure that contains files 4-6 (under the keys
  `tgm`,`tfs`, and `genes`) respectively as well as a `delim` key that describes
  what delimiter was used to separate out TFs in the case where there are
  multiple TFs in the same family.

- **events_to_genes_with_motifsregression_results.tsv**: Results from linear
  regression.

- **events_to_genes_with_motifsregression_results_FOREST_INPUT.tsv**: Only those
  regression results that fall under the p-value or q-value significance
  threshold provided in the configuration file, e.g. p=0.05, are included.
  This file can be used as input to Forest, and the prizes are -log2(pval)
  or -log2(qval).

- **regression_plots**: An optional subdirectory that contains plots visualizing
  the transcription factor linear regression tests.

Running forest.py
-----------------

Forest **requires** the compiled msgsteiner package.

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
                        for w, b, D. May contain values for optional
                        parameters mu, garnetBeta, noise, r, g. Default =
                        "./conf.txt"
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
  --musquared           Flag to add negative prizes to hub nodes proportional
                        to their degree^2, rather than degree. Must specify a
                        positive mu in conf file.
  --excludeTerms        Flag to exclude terminals when calculating negative
                        prizes. Use if you want terminals to keep exact
                        assigned prize regardless of degree.
  --msgpath=MSGPATH     Full path to the message passing code. Default =
                        "<current directory>/msgsteiner"
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
                        "_noisyEdges_" added to their names. The noise level
                        can be controlled using the configuration file.
                        Default = 0
  --shuffledPrizes=SHUFFLENUM
                        An integer specifying how many times you would like to
                        shuffle around the given prizes and re-run the
                        algorithm. Results of these runs will be merged
                        together and written in files with the word
                        "_shuffledPrizes_" added to their names. Default = 0
  --randomTerminals=TERMNUM
                        An integer specifying how many times you would like to
                        apply your given prizes to random nodes in the
                        interactome (with a similar degree distribution) and
                        re-run the algorithm. Results of these runs will be
                        merged together and written in files with the word
                        "_randomTerminals_" added to their names. Default = 0
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

### Forest input files and parameters

#### Required inputs

The first two options (`-p` and `-e`) are required. You should record your
terminal nodes and prize values in a text file. The file
`example/a549/Tgfb_phos.txt` is an example of what this file should look like.
You should record your interactome and edge weights in a text file with 3 or 4
columns. The file `data/iref_mitab_miscore_2013_08_12_interactome.txt` is a
human interactome example (this interactome comes from iRefIndex v13, scored and
formatted for our code).

A sample configuration file, `a549/tgfb_forest.cfg` is supplied. The user can
change the values included in this file or can supply their own
similarly formatted file. Unlike Garnet, the Forest configuration file name must
be preceded with `-c` or `--conf=`.
If the `-c` argument is not included in the command line
the program will attempt to read the default `conf.txt`. The parameters `w`, `b`, and `D`
must be set in the configuration file. Optional parameters `mu`, `garnetBeta`, `noise`,
`g`, and `r` may also be included.

```
w = float, controls the number of trees
b = float, controls the trade-off between including more
    terminals and using less reliable edges
D = int, controls the maximum path-length from v0 to terminal nodes
mu = float, controls the degree-based negative prizes (defualt 0.0)
garnetBeta = float, scales the garnet output prizes relative to the
             provided protein prizes (default 0.01)
noise = float, controls the standard deviation of the Gaussian edge
        noise when the --noisyEdges option is used (default 0.333)
g = float, msgsteiner parameter that affects the convergence of the
    solution and runtime (default 0.001)
r = float, msgsteiner parameter that adds random noise to edges,
    which is rarely needed because the Forest --noisyEdges option
    is recommended instead (default 0)
processes = int, number of processes to spawn when doing randomization runs
            (default to number of processors on your computer)

```

For more details about the parameters, see our publication.


#### Optional inputs

The rest of the command line options are optional.

If you have run the garnet module to create scores for transcription factors,
you can include that output file with the `--garnet` option and use `garnetBeta` in the
configuration file to scale the garnet scores.

The `--dummyMode` option will change which nodes in the terminal are connected
to the dummy node in the interactome. We provide an example of this using
`a549/Tgfb_interactors.txt`. For an explanation of the dummy node, see
publication.

The `--musquared` option will apply negative prizes to nodes based on their
squared degree, as opposed to linear degree. This is helpful if the default
mu behavior is not strict enough to eliminate irrelevant hub nodes from your
network.

If the file `msgsteiner` is not in the same directory as
forest.py, the path needs to be specified using the `--msgpath` option, e.g.,
'--msgpath /home/msgsteiner-1.3/msgsteiner'.

If you would like the output files to be stored in a directory other than the
one you are running the code from, you can specify this directory with the
`--outpath` option. The names of the output files will all start with the word
`result` unless you specify another word or phrase, such as an identifying label
for this experiment or run, with the `--outlabel option`. The `--cyto30` and
`--cyto28` tags can be used to specify which version of Cytoscape you would like
the output files to be compatiable with.

We include three options, `--noisyEdges`, `--shuffledPrizes`, and
`--randomTerminals` to determine how robust your results are by comparing them
to results with slightly altered input values. To use these options, supply a
number for either parameter greater than 0. If the number you give is more than
1, it will alter values and run the program that number of times and merge the
results together. The program will add Gaussian noise to the edge values you
gave in the `-e` option, or shuffle the prizes around all the network proteins
in the `-p` option, or assign the prizes to network proteins with similar
degrees as your original terminals, according to which option you use. In
`--noisyEdges`, Gaussian noise with mean 0 and standard deviation specified by
the parameter `noise` in the configuration file (default 0.333) will be added
to the edge scores. The results from these runs will be stored in seperate files
from the results of the run with the original prize or edge values, and both
will be outputted by the program to the same directory.

The knockout option can be used if you would like to simulate a knockout
experiment by removing a node from your interactome. Specify your knockout
proteins in a list, i.e. ['TP53'] or ['TP53', 'EGFR'].

The `-k` and `--cv` options can be used if you would like to run k-fold cross
validation. This will partition the proteins with prizes into k equal
subsamples. It will run msgsteiner k times, leaving one subsample of prizes out
each time. The `--cv-reps` option can be used if you would like to run k-fold
cross validation multiple times, each time with a different random partitioning
of terminals. If you do not supply `--cv-reps` but do provide a k, cross
validation will be run once. Each time it is run, a file called
`<outputlabel>_cvResults_<rep>.txt` will be created. For each of the k
iterations, it will display the number of terminals held out of the prizes
dictionary, the number of those that were recovered in the optimal network as
Steiner nodes, and the total number of Steiner nodes in the optimal network.

The `-s` option will supply a seed option to the pseudo-random number generators
used in noisyPrizes, shuffledPrizes, randomTerminals, and the optimization in
msgsteiner itself. If you want to reproduce exact results, you should supply the
same seed every time. If you do not supply your own seed, system time is used a
seed.

###Running forest

Once you submit your command to the command line the program will run. It will
display messages as it completes, letting you know where in the process you are.
If there is a warning or an error it will be displayed on the command line. If
the run completes successfully, several files will be created. These files can
be imported into Cytoscape v.3.0 to view the results of the run.  These files
will be named first with the outputlabel that you provided (or `result` by
default), and then with a phrase identifying which file type it is.

### Forest output

- **info.txt** contains information about the algorithm run, including any error
  messages if there were any during the run.

- **optimalForest.sif** contains the optimal network output of the
  message-passing algorithm (without the dummy node). It is a Simple Interaction
Format file. To see the network, open Cytoscape, and click on File > Import >
Network > File..., and then select this file to open. Click OK.

- **augmentedForest.sif** is the same thing, only it includes all the edges in
  the interactome that exist between nodes in the optimal Forest, even those
edges not chosen by the algorithm. Betweenness centrality for all nodes was
calculated with this network.

- **dummyForest.sif** is the same as optimalForest.sif, only it includes the
dummy node and all edges connecting to it.

- **edgeattributes.tsv** is a tab-seperated value file containing information
  for each edge in the network, such as the weight in the interactome, and the
fraction of optimal networks this edge was contained in. To
import this information into Cytoscape, first import the network .sif file you
would like to view, and then click on File > Import > Table > File..., and
select this file. Specify that this file contains edge attributes, rather than
node attributes, and that the first row of the file should be interpreted as
column labels. Click OK.

- **nodeattributes.tsv** is a tab-seperated value file containing information
  for each node in the network, such as the prize you assigned to it and
betweenness centrality in the augmented network. To import this information into
Cytoscape, first import the network .sif file you would like to view, and then
click on File > Import > Table > File..., and select this file. Specify that
this file contains node attributes, rather than edge attributes, and that the
first row of the file should be interpreted as column labels. Click OK.

When the network and the attributes are imported into Cytoscape, you can alter
the appearance of the network as you usually would using VizMapper.

Testing
-----------------
See the `tests` directory for instructions on testing Omics Integrator.
