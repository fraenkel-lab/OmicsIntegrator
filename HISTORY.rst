.. :changelog:

History
-------

0.3.1 (2017-02-02)
---------------------
* Garnet integration test case using FOXA1 ChIP-seq
* Correct formatting problems in example data wgEncodeRegDnaseClusteredV2.fasta.gz
* Fix bug when reading Garnet results in Forest
* Use relative paths for motif gifs in Garnet HTML output file
* Update sort order in Garnet output table
* Garnet uses balanced binary search tree to speed up the search for genes within peaks
* Change default values of arguments in Garnet's map_peaks_to_known_genes.py

0.3.0 (2016-09-03)
---------------------
* Remove conflicts between Forest edge noise and garnetBeta parameters
* Correctly query for whether prize nodes are in a directed network
* Use correct p-value or q-value to filter Garnet output
* Switch to msgsteiner-1.3 to fix bugs with directed edges and suboptimal forests
* BSD license
* Forest unit tests to assess the effects of w, beta, mu, and D
* Continuous integration testing with Travis CI
* Option to save plots from Garnet motif regression
* Write Forest objective function value of optimal forest to info file

0.2.0 (2016-04-13)
---------------------
* Add a Forest integration test
* Correct nondeterministic msgsteiner behavior
* Add additional example data sets and scripts
* Switch to Creative Commons Attribution-NonCommercial 4.0 International Public License
* Rename Forest configuration file parameter n to garnetBeta
* Add Forest configuration file parameter processes
* Correct error in negative prize score calculation
* Improve Garnet output file format
* Add Forest option to exclude terminals

0.1.0 (2015-04-01)
---------------------
* First release on PyPI.
