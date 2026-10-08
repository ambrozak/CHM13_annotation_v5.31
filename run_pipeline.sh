#!/bin/bash

# Create mapping from GRCh38 MANE to CHM13:

lifton -g ./MANEv1.5.gff -chroms ./chrom_mapping.txt -copies -sc 0.9 ./chm13v2.0.fa ./hg38.p14.fa -o MANE_on_chm13.gff

# Give each pseudogene in CHM13v5.3 a transcript using gffread's --gene2exon flag:

gffread --gene2exon ./chm13v2.0.fa > ./chm13v2.0_gene2exon.gff

# Compare annotations using gffcompare:

gffcompare -r ./chm13v2.0_gene2exon.gff ./MANE_on_chm13.gff

# Run the comparison script:

python compare_proteins.py

# Run the annotation-merging script:

python update_annotation.py