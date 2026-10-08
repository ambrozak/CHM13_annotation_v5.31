#!/bin/bash

# Download CHM13 genome

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

gunzip chm13v2.0.fa.gz

# Download CHM13v5.3 Annotation

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz

gunzip chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz

mv chm13v2.0_RefSeq_Liftoff_v5.3.gff chm13v5.3.gff

# Download GRCh38 genome

wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/p14/hg38.p14.fa.gz

gunzip hg38.p14.fa.gz


# Download MANEv1.5 Annotation

wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.refseq_genomic.gff.gz

gunzip MANE.GRCh38.v1.5.refseq_genomic.gff.gz 

mv MANE.GRCh38.v1.5.refseq_genomic.gff MANEv1.5.gff