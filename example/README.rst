===============================
garnet+forest example files
===============================

This directory includes sample datasets that can be used to run your own data
analysis or test out the scripts with existing data.
           

Files
-----
*a549/
    Here we have data from lung cancer-derived cell lines that we use to evaluate the role of TgfB-stimulation on gene expression and phospho-proteomic changes.
    --tgfb_garnet.cfg: Configuration file required to use garnet to identify transcription factors from epigenetic data and expression data.
    --Tgfb_phos.txt: Changes in phospho-protein levels between TgfB-stimulated cells and untreated cells.
    --Tgfb_exp.txt:  Changes in gene expression levels between TgfB-stimulated cells and untreated cells.
    --tgfb_interactors.txt: Proteins that interact with TgfB, used to run forest algorithm.
    --wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.bed: BED-formatted file of DNase-I hypersensitive regions in a549 cell lines.
    --wgEncodeUWDukeDnaseA549.fdr01peaks.hg19.fasta.gz: FASTA-formated file of the same regions.

*mcf7/
    Here we have data from a breast cancer-derived cell line that we can use to construct TF binding scores that can be used with any breast-cancer related gene expression dataset of choice.
    --mcf7_garnet.cfg
    --wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.bed
    --wgEncodeUWDukeDnaseMCF7.fdr01peaks.hg19.fasta.gz

*dnaseClus/
    Here we have clustered DNase-I hypersensitive regions across the human genome. This data can be used in cases where the user has no specific epigenetic data. 
    --dnaseClus_garnet.cfg:
    --wgEncodeRegDnaseClusteredV2.bed:
    --wgEncodeRegDnaseClusteredV2.fasta.gz:

*murineFib/
    --murineFib_garnet.cfg
    --wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.narrowPeak
    --wgEncodeUwDnaseFibroblastC57bl6MAdult8wksPk_Rep1AND2.fasta.gz

