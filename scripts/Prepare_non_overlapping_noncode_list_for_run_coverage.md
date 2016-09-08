
---
title: "Prepare a list of non-overlapping noncodes and save it in bed format"
author: "Priscila Darakjian"
output: html_document
---

# Prepare a list of non-overlapping noncodes and save it to a bed file
#### author: Priscila Darakjian


### Required Libraries

```r
library(GenomicFeatures)
library(Rsamtools)
```

### Read the noncode gtf annotation file and extract annotation for exons only

```r
setwd("/lawrencedata/ongoing_analyses/RNASeq016/RNASeq016_noncode/NonCode_Analysis_2016")
gtf<-read.delim("/lawrencedata/ongoing_analyses/RNASeq016/RNASeq016_noncode/NonCode_Analysis_2016/data/NONCODE2016_mouse_mm10_lncRNA.gtf", sep="\t", header=FALSE, stringsAsFactors=FALSE)
gtf <- gtf[gtf$V3 == "exon",]
```

### Extract the gene symbol id field (1) from column 9 of gtf and remove the string "gene_id " from it

```r
sym.gene.name.gtf<-sapply(strsplit(gtf$V9,";"),function(x) x[1])
sym.gene.name.gtf<-sub("gene_id ","",sym.gene.name.gtf)
```

### Create a separate column for noncode ids and change where "strand" field (7) is "." to "*"

```r
gtf$gene.sym<-sym.gene.name.gtf
gtf$V7[gtf$V7 == "."] <- "*"
```

### Generate a GRanges object from gtf above for working with chr ranges, and find overlaps between exons

```r
gtf.ranges<-GRanges(seqnames=Rle(gtf$V1),ranges=IRanges(start=gtf$V4,end=gtf$V5),strand=gtf$V7,gene.names=gtf$gene.sym)
gtf.range.ovls<-findOverlaps(gtf.ranges,ignoreSelf=T,ignoreRedundant=T,type="any")
gtf.range.ovls.mat <- as.matrix(gtf.range.ovls)
gtf.range.ovls.dta <- data.frame(as.data.frame(gtf.ranges)[gtf.range.ovls.mat[,1],], as.data.frame(gtf.ranges)[gtf.range.ovls.mat[,2],])
```

### Verify if there are overlaps between exons of different genes and prepare a list of the ranges in which exons are overlapping between genes

```r
gtf.ovl.genes.inds <- which(gtf.range.ovls.dta$gene.names != gtf.range.ovls.dta$gene.names.1)
gtf.rm.inds <- as.integer(gtf.range.ovls.mat[gtf.ovl.genes.inds,])
gtf.ranges.gene.ovl.df <- as.data.frame(gtf.ranges[unique(gtf.rm.inds)])
```

### Convert data frame to bed format so that we can run bedtools (intersectBed and subtractBed)

```r
gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[c(1,2,3,6,4,5)]
```
***
##### *This wasn't run --> Remove lines with non-standard chromosome names*
##### *gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[- grep("random", gtf.ranges.gene.ovl.df$seqnames),]*
##### *gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[- grep("chrUn_", gtf.ranges.gene.ovl.df$seqnames),]*
***

### This is the original list of exons

```r
gtf.ranges.df.unique<-unique(as.data.frame(gtf.ranges))
gtf.ranges.df.unique<-gtf.ranges.df.unique[c(1,2,3,6,4,5)]
```

### Generate bed files for exons overlapping between genes and for the original list of exons 

```r
write.table(gtf.ranges.gene.ovl.df,"overlapping_exons_btwn_genes.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(gtf.ranges.df.unique,"exons.bed",quote=F,col.names=F,row.names=F,sep="\t")
```

### This is run in Linux to remove from original list just the portions of exons that overlap between genes
#### sort -k1,1 -k2,2n overlapping_exons_btwn_genes.bed > overlapping_exons_btwn_genes_sorted.bed
#### sort -k1,1 -k2,2n exons.bed > exons_sorted.bed
#### /lawrencedata/darakjia/bedtools-2.17.0/bin/intersectBed -a exons_sorted.bed -b overlapping_exons_btwn_genes_sorted.bed > intersect.bed
#### /lawrencedata/darakjia/bedtools-2.17.0/bin/subtractBed -s -a exons_sorted.bed -b intersect.bed > non_overlapping_exons.bed

#*****************************************************************************
# IN THE CASE OF MF5 (CYNO) THERE WERE NO OVERLAPPING EXONS BETWEEN DIFFERENT GENES
# SO:
# cat exons_sorted.bed  > non_overlapping_exons.bed 
#*****************************************************************************
# Get back the file with no gene overlaps
no_gene_ovl_exons <- read.table("non_overlapping_exons.bed")
# Convert it to genomic ranges
no.gene.ovl.exons.ranges<-GRanges(seqnames=Rle(no_gene_ovl_exons$V1),ranges=IRanges(start=no_gene_ovl_exons$V2,end=no_gene_ovl_exons$V3),strand=no_gene_ovl_exons$V6,gene.names=no_gene_ovl_exons$V4)

# Merge exon isoforms
no.gene.ovl.exons.ranges.red<-reduce(split(no.gene.ovl.exons.ranges,elementMetadata(no.gene.ovl.exons.ranges)$gene.names))

# Convert above objects into data frame
no.gene.ovl.exons.ranges.red.df<-as.data.frame(no.gene.ovl.exons.ranges.red)

# Re-arrange data frame in bed format and write object to file
no.gene.ovl.exons.ranges.red.df<-no.gene.ovl.exons.ranges.red.df[c(3,4,5,2,6,7)]
no.gene.ovl.exons.ranges.red.df$width<-0

write.table(no.gene.ovl.exons.ranges.red.df,"MF5_reduced_non_overlapping_exons.txt", sep="\t",quote=F,col.names=F,row.names=F)

#Run coverage with bedtools
#sort -k6 -k1,1 -k2,2n UNMC_reduced_non_overlapping_exons.txt > UNMC_reduced_non_overlapping_exons_sorted.bed
#sed "s/^/chr/g" UNMC_reduced_non_overlapping_exons_sorted.bed > UNMC_reduced_non_overlapping_exons_sorted_chradd.bed
#nohup ./run_bedtools_coverage.sh > run_bedtools_coverage.log 2>&1 </dev/null &

#NOTE: I HARDCODED THE NAMES OF THE SAMPLES IN run_bedtools_coverage.sh
#      to make sure we know the order of the columns in the read coverage result file.
#      I ordered the sample names by doing:
#      tr ',\n' ' ' < bam_files.txt > run_bedtools_coverage2.sh 
#      and then completing the sh file with the code around those names.
