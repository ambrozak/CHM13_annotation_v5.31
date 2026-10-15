# T2T-CHM13 Gene Annotation v5.31

## Motivation

This update incorporates newly added genes from MANE v1.5 into the T2T-CHM13 annotation. We mapped MANE v1.5 from GRCh38 to CHM13 using LiftOn, which refines LiftOff-based mappings and can identify additional copies of genes. This yielded both improved mappings for genes already annotated in CHM13 and mappings in novel loci.

## Results

The v5.31 annotation update adds 16 novel transcripts to existing genes, 8 improvements to existing transcript mappings, 17 novel gene loci (10 of which were extra copies of existing genes), 2 new pseudogene annotations, and replaces 2 pseudogene annotations with genes. 

A count of gene biotypes before and after can be found below:


| Gene biotype | v5.3 | v5.31 |
|----------|----------|----------|
| Protein coding   | 20,064    | 20,083 (+19)    |
| LncRNA    | 18,385     | 18,385 (+0)    |
| Pseudogene    | 17,466  | 17,466 (+2, -2)   |
| Other	| 5,427 | 5,427 (+0)	|
| **Total** | 61,342 | 61,361 (+19) |

Note these totals are counted using the gene_biotype field, differing slightly from the [previous release](https://github.com/oldvalley49/chm13_annotation_v5.3/tree/master).

## Methods

### Lift-Over and Mapping

We first ran LiftOn, mapping over MANE v1.5 genes from GRCh38 onto CHM13 using the following command:

```bash
lifton -g ./MANEv1.5.gff -chroms ./chrom_mapping.txt -copies -sc 0.9 ./chm13v2.0.fa ./hg38.p14.fa -o MANE_on_chm13.gff
```

This lifts chromosome-by-chromosome, and searches for additional copies of genes that share at least 90% identity with the first copy.

We altered the CHM13 v5.3 annotation to give each pseudogene a single-exon transcript, then used gffcompare to map each MANE transcript to its closest CHM13 v5.3 transcript:

```bash
gffread --gene2exon ./chm13v2.0.fa > ./chm13v2.0_gene2exon.gff

gffcompare -r ./chm13v2.0_gene2exon.gff ./MANE_on_chm13.gff
```

### Filtering Out Low Quality Mappings

After running gffcompare, we filtered down the MANE gene mappings to a set that are likely to be 1) correctly mapped, and 2) an improvement over the existing CHM13 annotations in the same locus, if there are any.

First, we filtered out MANE gene mappings that were likely to be incorrectly mapped based on gene order. We consider each gene's immediate neighbors both before and after mapping. We retain a mapping only if at least one of its original neighbors remains adjcent to it after mapping. For example, take this gene ordering:

```
Before mapping:
 ... - A - B - C - D - E - F - ...

After mapping:
 ... - A - E - C - D - B - F - ...
```

Here, E and B may have switched positions due to high similarity. We automatically filter out both because neither of their neighbors stay the same after mapping. We still consider C and D because at least one of their neighbors is conserved.

Second, to filter out MANE mappings that are unlikely to improve the existing annotation, we use the following criteria:

```
MANE transcript A is paired with a transcript on an existing CHM13 gene G by gffcompare.

If G contains a transcript with the same transcript ID as A, select that transcript as B.
Otherwise, select the transcript on G assigned to A by gffcompare as B.

If A does not share a db_xref attribute with B:
    Filter out A

If both A and B are protein coding:
    If B's CDS has an equal or higher protein identity to A's original CDS on GRCh38 than A's CDS:
        Filter out A
```

After this process, we're left with a set of MANE gene mappings that are likely to be improvements to existing coding annotations (33), in novel loci (17), or overlap existing pseudogenes (143).

### Replacing or Adding To Existing Annotations

We then decided whether to include or exclude the new mappings using the following logic:

```
Given MANE transcript A on CHM13


If A is matched with existing coding transcript B:
    If B's CDS has at least 80% protein identity with A's original MANE protein:
        If B has the same transcript ID as A:
            Automatically replace B with A
        If B has a different transcript ID:
            Automatically add A
    Else:
        Manually review A


If A is in a novel locus:
    If A's CDS has at least 95% protein identity with A's original MANE protein:
        If A has no missing start/stop codon or frameshift mutations:
            Automatically add A
    Else:
        Manually review A


If A overlaps with an existing pseudogene P:
    If A's CDS has 100% protein identity with A's original MANE protein:
        If LiftOn did not find extra copies of A on CHM13 which overlap with an existing functional annotation of A:
            Automatically replace P with A
    Else:
        Filter out A
```

23 mappings were manually reviewed, final decisions and notes on each can be found in [manual_review_decisions.csv](./manual_review_decisions.csv)

## Reproducing Results

To run these scripts, install the following:

Tools:
 - LiftOn
 - gffcompare
 - gffread

Python Libraries:
 - pandas
 - pyfaidx
 - biopython
 - tqdm

 Afterwards, simply run:

 ```bash
 bash download_data.sh

 bash run_pipeline.sh
 ```

 Note that generating the new annotation relies on the decisions in [manual_review_decisions.csv](./manual_review_decisions.csv).

 ### Output

 The updated annotation file is created as ```chm13v5.31_MANE_updated.gff```. 
 A list of all changes made to the annotation is in ```annotation_changes.tsv```. 
 Lastly, a record of what was done with the novel MANE v1.5 genes is listed in ```gene_trace.tsv```.

 If you wish to track a different set of genes from MANE, change the list of dbxref gene IDs in ```gene_ids_new_in_MANE1.5.txt``` to those of your choosing.