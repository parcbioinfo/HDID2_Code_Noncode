
---
title: "Prepare a list of non-overlapping code and noncode genes and gap intervals and save it in bed format"
author: "Priscila Darakjian"
output: html_document
---
# Script: Prepare_non_overlapping_code-and-noncode_list_for_running_coverage.Rmd
## Description: Prepare a list of non-overlapping code and noncode genes and save it to a bed file in order to obtain read coverage for those genomic intervals. Also prepare a file containing the gap intervals (complement of the code-noncode bed file)
#### Author: Priscila Darakjian - OHSU - Research - Behavioral Neuroscience
#### Date: 09/22/2016

### Required Libraries

```r
library(GenomicFeatures)
library(Rsamtools)

load("data/Combined_annot_clean.RData")

setwd("/lawrencedata/ongoing_analyses/RNASeq016/RNASeq016_noncode/NonCode_Analysis_2016")
```
### Convert annotation into genomic ranges

```r
Comb_rgs <- GRanges(seqnames = Rle(Combined_annot_clean$seqnames), 
    ranges = IRanges(start = Combined_annot_clean$start, end = Combined_annot_clean$end), 
    strand = Combined_annot_clean$strand, gene.names = Combined_annot_clean$gene.names, 
    biotype = Combined_annot_clean$biotype)
```

### Find overlapping intervals. Ultimately we want to remove the sections of those intervals which overlap between different genes.

```r
Comb_rgs_ovls <- findOverlaps(Comb_rgs, drop.self = T, drop.redundant = T, 
    type = "any")
Comb_rgs_ovls.mat <- as.matrix(Comb_rgs_ovls)
Comb_rgs_ovls.dta <- data.frame(as.data.frame(Comb_rgs)[Comb_rgs_ovls.mat[, 
    1], ], as.data.frame(Comb_rgs)[Comb_rgs_ovls.mat[, 2], ])
```

### Verify if there are overlaps between exons of different genes and prepare a list of the ranges in which exons are overlapping between genes

```r
Comb_rgs_ovls.genes.inds <- which(Comb_rgs_ovls.dta$gene.names != 
    Comb_rgs_ovls.dta$gene.names.1)
Comb_rgs_ovls.rm.inds <- as.integer(Comb_rgs_ovls.mat[Comb_rgs_ovls.genes.inds, 
    ])
Comb_rgs.gene.ovl.df <- as.data.frame(Comb_rgs[unique(Comb_rgs_ovls.rm.inds)])
```

### Convert data frame to bed format so that we can run bedtools (intersectBed and subtractBed)

```r
descriptor <- as.data.frame(paste(Comb_rgs.gene.ovl.df$gene.names, 
    ":", Comb_rgs.gene.ovl.df$biotype))
Comb_rgs.gene.ovl.df = cbind(Comb_rgs.gene.ovl.df, descriptor)
Comb_rgs.gene.ovl.bed <- Comb_rgs.gene.ovl.df[c(1, 2, 3, 8, 4, 
    5)]
```

### This is the original list of exons

```r
Comb_rgs.df <- as.data.frame(Comb_rgs)
descriptor <- as.data.frame(paste(Comb_rgs.df$gene.names, ":", 
    Comb_rgs.df$biotype, sep = ""))
Comb_rgs.df = cbind(Comb_rgs.df, descriptor)
Comb_rgs.bed <- Comb_rgs.df[c(1, 2, 3, 8, 4, 5)]
```

### Generate bed files for exons overlapping between genes and for the original list of exons 

```r
write.table(Comb_rgs.gene.ovl.bed, "data/overlapping_exons_btwn_genes.bed", 
    quote = F, col.names = F, row.names = F, sep = "\t")
write.table(Comb_rgs.bed, "data/exons.bed", quote = F, col.names = F, 
    row.names = F, sep = "\t")
```

### This is run in Linux to remove from original list just the portions of exons that overlap between genes
#### /lawrencedata/darakjia/bedtools-2.17.0/bin/intersectBed -s -a exons_sorted.bed -b overlapping_exons_btwn_genes_sorted.bed > intersect.bed
#### /lawrencedata/darakjia/bedtools-2.17.0/bin/subtractBed -s -a exons_sorted.bed -b intersect.bed > non_overlapping_exons.bed

### Get back the file with no gene overlaps and convert it to genomic ranges

```r
no_gene_ovl_exons <- read.table("data/non_overlapping_exons.bed")
no.gene.ovl.exons.ranges <- GRanges(seqnames = Rle(no_gene_ovl_exons$V1), 
    ranges = IRanges(start = no_gene_ovl_exons$V2, end = no_gene_ovl_exons$V3), 
    strand = no_gene_ovl_exons$V8, gene.names = no_gene_ovl_exons$V4, 
    biotype = no_gene_ovl_exons$V6)
```

### Merge exon isoforms

```r
no.gene.ovl.exons.ranges.red <- reduce(split(no.gene.ovl.exons.ranges, 
    elementMetadata(no.gene.ovl.exons.ranges)$gene.names))

# Convert above reduced ranges into a data frame and add the
# biotype information to it
no.gene.ovl.exons.ranges.red.df <- as.data.frame(no.gene.ovl.exons.ranges.red)
no.gene.ovl.exons.ranges.red.df.ordered <- no.gene.ovl.exons.ranges.red.df[order(no.gene.ovl.exons.ranges.red.df$strand, 
    no.gene.ovl.exons.ranges.red.df$seqnames, no.gene.ovl.exons.ranges.red.df$start), 
    ]
comb_mdata <- unique(Combined_annot_clean[, 6:7])
no.gene.ovl.exons.ranges.red.df.ordered$gene.names = no.gene.ovl.exons.ranges.red.df.ordered$group_name
no.gene.ovl.exons.ranges.red.df.ordered <- no.gene.ovl.exons.ranges.red.df.ordered[, 
    3:8]
no.gene.ovl.exons.ranges.red.df.ordered$biotype = comb_mdata[match(no.gene.ovl.exons.ranges.red.df.ordered$gene.names, 
    comb_mdata$gene.names), "biotype"]
```
### Re-arrange data frame in bed format  

```r
descriptor <- as.data.frame(paste(no.gene.ovl.exons.ranges.red.df.ordered$gene.names, 
    ":", no.gene.ovl.exons.ranges.red.df.ordered$biotype, sep = ""))
no.gene.ovl.exons.ranges.red.df.ordered = cbind(no.gene.ovl.exons.ranges.red.df.ordered, 
    descriptor)
no.gene.ovl.exons.ranges.red.df.bed <- no.gene.ovl.exons.ranges.red.df.ordered[c(1, 
    2, 3, 8, 4, 5)]
```
### Replace "width" with zeros, and add "chr" to the beginning of each "seqname" (chromosome) to match those in the bam files and write object to file 

```r
no.gene.ovl.exons.ranges.red.df.bed$width <- 0
no.gene.ovl.exons.ranges.red.df.bed$seqnames <- paste("chr", 
    no.gene.ovl.exons.ranges.red.df.bed$seqnames, sep = "")
write.table(no.gene.ovl.exons.ranges.red.df2.bed, "data/mm10_noncode_reduced_non_overlapping_exons.bed", 
    sep = "\t", quote = F, col.names = F, row.names = F)
```

### Generate temporary separate strand files in order to obtain the "complement" intervals (gaps) through bedtools for each strand.

```r
pos_strand <- no.gene.ovl.exons.ranges.red.df.bed[no.gene.ovl.exons.ranges.red.df.bed$strand == 
    "+", ]
neg_strand <- no.gene.ovl.exons.ranges.red.df.bed[no.gene.ovl.exons.ranges.red.df.bed$strand == 
    "-", ]
write.table(pos_strand, "data/mm10_noncode_reduced_non_overlapping_exons_pos_strand.bed", 
    sep = "\t", quote = F, col.names = F, row.names = F)
write.table(neg_strand, "data/mm10_noncode_reduced_non_overlapping_exons_neg_strand.bed", 
    sep = "\t", quote = F, col.names = F, row.names = F)
```
### Run the following commands in bedtools (NOTE: for some chromosomes I had to increase the length value since some noncodes exceed the chromosome length for some reason)
#### bedtools complement -i mm10_noncode_reduced_non_overlapping_exons_neg_strand.bed -g mm10.genome > complement_neg.bed
#### bedtools complement -i mm10_noncode_reduced_non_overlapping_exons_pos_strand.bed -g mm10.genome > complement_pos.bed

### Reload the temporary files to merge them and add strand information

```r
compl_pos <- read.table("data/complement_pos.bed")
compl_pos$feature <- "gap"
compl_pos$score <- "0"
compl_pos$strand <- "+"
compl_neg <- read.table("data/complement_neg.bed")
compl_neg$feature <- "gap"
compl_neg$score <- "0"
compl_neg$strand <- "-"
```
### This is the file that will be used to generate the coverage for the gaps

```r
gaps_bed <- rbind(compl_pos, compl_neg)
write.table(gaps_bed, "mm10_gaps.bed", sep = "\t", quote = F, 
    col.names = F, row.names = F)
```

### Run coverage with bedtools for both files (combined code/noncode and gaps)
#### nohup ./run_bedtools_coverage.sh > run_bedtools_coverage.log 2>&1 </dev/null &
#### nohup ./run_bedtools_code-gaps_coverage.sh > run_bedtools_coverage.log 2>&1 </dev/null &
##### NOTE: I HAVE HARDCODED THE NAMES OF THE SAMPLES IN run_bedtools_coverage.sh (rather than using a wildcard) to make sure we know the order of the columns in the read coverage result file. I have created a file with the sample names (ls *.bam > bam_files.txt) and added those to the script that will run bedtoos mutibamcov by issuing:
###### tr ',\n' ' ' < bam_files.txt > run_bedtools_coverage2.sh 
##### and then I completed the shell script with the bedtools command around those names.
