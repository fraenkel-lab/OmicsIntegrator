===============================
Omics Integrator Examples
===============================

This directory includes sample datasets **and test scripts** that can be used to run your own data
analysis or test out the scripts with existing data. We divided the examples based on tissue
of origin.

Files in the `a549` directory will enable you to run Garnet + Forest on your system.  `test-tgfb-data.py` is a script that will run both Garnet + Forest, or show you how to just run Forest on a single dataset.

If you want to test out Omics Integrator without providing epigenetic data, then you can run `build-dnase-matrix.py` in the `dnaseClust` directory (for human) or `build-murine-fib-matrix.py` in the `murineFib` directory (for mouse).


a549/
----
Here we have data from lung cancer-derived cell lines that we use to evaluate the role of TgfB-stimulation on gene expression and phospho-proteomic changes. Running `test-tgfb-data.py` at the prompt will first run garnet.py, then run Forest with the result.

- A549_FOXA1_broadPeak.bed: [FOXA1 ChIP-seq data from ENCODE](https://www.encodeproject.org/experiments/ENCSR000BPX/)
- A549_FOXA1_broadPeak.fasta: [FOXA1 ChIP-seq data from ENCODE](https://www.encodeproject.org/experiments/ENCSR000BPX/)
- test-tgfb-data.py: A test script designed to show how to use Garnet + Forest to interpret real data.
- tgfb_forest.cfg: Configuration file require to use Forest to find interactions between altered proteins during TgfB stimulation.
- tgfb_foxa1_garnet.cfg: Configuration file required to test Garnet with FOXA1 ChIP-seq data.
- tgfb_garnet.cfg: Configuration file required to use Garnet to identify transcription factors from epigenetic data and expression data.
- Tgfb_phos.txt: Changes in phospho-protein levels between TgfB-stimulated cells and untreated cells.
- Tgfb_exp.txt:  Changes in gene expression levels between TgfB-stimulated cells and untreated cells.
- tgfb_interactors.txt: Proteins that interact with TgfB, used to run Forest.
- wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed: BED-formatted file of DNase-I hypersensitive regions in a549 cell lines.
- wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.fasta.gz: FASTA-formated file of the same regions.

mcf7/
----
Here we have data from a breast cancer-derived cell line that we can use to construct TF binding scores that can be used with any breast-cancer related gene expression dataset of choice.

- build-mcf7-matrix.py: This script shows how to use Garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- mcf7_garnet.cfg: Configuration file used.
- wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.bed: DNase-I hypersensitive regions in MCF7 cell lines.
- wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.fasta.gz: FASTA sequences derived from the above BED file.

dnaseClus/
---------
Here we have clustered DNase-I hypersensitive regions across the human genome. This data can be used in cases where the user has no specific epigenetic data.

- build-dnase-matrix.py: This script shows how to use Garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- dnaseClus_garnet.cfg: Configuration file used.
- wgEncodeRegDnaseClusteredV2.bed: DNase-I hypersensitive regions clustered across all human cell lines in ENCODE.
- wgEncodeRegDnaseClusteredV2.fasta.gz:FASTA sequences derived from the above BED file.

murineFib/
----------
Here we have DNase-I hypersensitive data from a murine fibroblast cell line. This data can be used
for mouse expression data.

- build-murin-fib-matrix.py: This script shows how to use Garnet to predict transcription factor binding from epigenetic data (without regression to select best transcription factors).
- murineFib_garnet.cfg: Configuration file used.
- wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.narrowPeak: DNase-I hypersensitive regions from murine fibroblasts.
- wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.fasta.gz: FASTA sequences derived from the above BED file.

CPDB/
---------
Here we showcase the ability of the forest functionality to recapitulate known
pathways from the Consensus Pathway Database better than the tree formulation of
Forest. This example also demonstrates one strategy for running Forest with
many different parameter combinations.

- CPDB_case_study.py: Script to run the example. You must replace the `msgpath` variable with the path to the msgsteiner executable on your system.
- pyr_rho_mrna_noUBC.terminal: Terminal file to use

GBM/
---------
Here we demonstrate in silico knockouts in Forest using phosphoproteomic data from glioblastoma U87 cell lines.

- GBM_case_study.py: Script to run Forest with the full network and with a simulated EGFR knockout. You must replace the `msgpath` variable with the path to the msgsteiner executable on your system.
- gbm_prize.txt: Terminal file to use
