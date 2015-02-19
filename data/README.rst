====================
garnet + forest data
====================

This directory includes basic data files required to run garnet and
forest. While all aspects are fully configurable, we provide the files
to facilitate using garnet and forest without additional downloads.

Files
-----
- **matrix_files/**
This directory contains the TAMO-formatted motif files required to
scan epigenetic regions
- **kgXrefFiles:**
These files map UCSC identifiers to human and mouse gene names.
  - ucsc_hg19_kgXref.txt: Annotations for hg19
  - ucsc_mm9_kgXref.txt: Annotations for mm9

- **KnownGene files:**
These files contain UCSC known genes and are required to map epigenetic regions genes so 
that predicted transcription factor binding sites can be associated with changes in gene
expression
  - ucsc_hg19_knownGenes.txt: Annotations for hg19
  - ucsc_mm9_knownGenes.txt: Annotations for mm9
- **Protein-protein interaction network:**
  We included a recent version of the human protein-protein interaction network to facilitate
running forest.  
  - iref_mitab_miscore_2013_08_12_interactome.txt: We collected the IRefWeb interactome using
the MIScore function.

