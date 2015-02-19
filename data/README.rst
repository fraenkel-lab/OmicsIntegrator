===============================
data required to run garnet and forest
===============================

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
  This includes ucsc_hg19_kgXref.txt and ucsc_mm9_kgXref.txt
- **KnownGene files:**
  These files contain UCSC known genes and are required to map epigenetic regions genes so 
  that predicted transcription factor binding sites can be associated with changes in gene
  expression
- **Protein-protein interaction network:**
  We included a recent version of the human protein-protein interaction network to facilitate
  running forest.  

