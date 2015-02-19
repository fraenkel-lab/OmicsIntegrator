===============================
garnet+forest examples
===============================

This directory includes sample datasets **and test scripts** that can be used to run your own data
analysis or test out the scripts with existing data. We divided the examples based on tissue
of origin. 
           

a549/
----
Here we have data from lung cancer-derived cell lines that we use to evaluate the role of TgfB-stimulation on gene expression and phospho-proteomic changes. Running `test-tgfb-data.py` at the prompt will first run garnet.py, then run forest with the result.

- test-tgfb-data.py: A test script designed to show how to use garnet + forest to interpret real data.
- tgfb_forest.cfg: Configuration file require to use forest to find interactions between altered proteins during TgfB stimulation.
- tgfb_garnet.cfg: Configuration file required to use garnet to identify transcription factors from epigenetic data and expression data.
- Tgfb_phos.txt: Changes in phospho-protein levels between TgfB-stimulated cells and untreated cells.
- Tgfb_exp.txt:  Changes in gene expression levels between TgfB-stimulated cells and untreated cells.
- tgfb_interactors.txt: Proteins that interact with TgfB, used to run forest algorithm.
- wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed: BED-formatted file of DNase-I hypersensitive regions in a549 cell lines.
- wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.fasta.gz: FASTA-formated file of the same regions.

mcf7/
----
Here we have data from a breast cancer-derived cell line that we can use to construct TF binding scores that can be used with any breast-cancer related gene expression dataset of choice.

- test-mcf7-data.py: This script shows how to use garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- mcf7_garnet.cfg
- wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.bed
- wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.fasta.gz

dnaseClus/
---------
Here we have clustered DNase-I hypersensitive regions across the human genome. This data can be used in cases where the user has no specific epigenetic data. 

- test-dnase-data.py: This script shows how to use garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- dnaseClus_garnet.cfg:
- wgEncodeRegDnaseClusteredV2.bed:
- wgEncodeRegDnaseClusteredV2.fasta.gz:

murineFib/
----------
Here we have DNase-I hypersensitive data from a murine fibroblast cell line. This data can be used
for mouse expressiond data.

- test-murine-data.py: This script shows how to use garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- murineFib_garnet.cfg
- wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.narrowPeak
- wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.fasta.gz

