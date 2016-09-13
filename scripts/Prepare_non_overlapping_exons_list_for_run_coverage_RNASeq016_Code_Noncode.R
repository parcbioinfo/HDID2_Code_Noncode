library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)

load("noncode.cln.RData")
Ensembl <-read.table("gtfs/Mus_musculus.NCBIM37.67_exons.gtf",sep="\t", header=F,stringsAsFactors=F)
# Uncomment below if noncode.RData is not available (and if file below is available)
# Noncode <-read.table("NONCODEv4_mouse_lncRNA_exons_nonredundant.txt",header=T,stringsAsFactors=F)
Noncode <- noncode.cln
#Noncode<-Noncode[- grep("MT_", Noncode$seqnames),]
Ensembl <- Ensembl[Ensembl$V3 == "exon",]

# Extract gene symbol field (4 and 1, respectively) from column 9 of gtfs
ens.gene.name<-sapply(strsplit(Ensembl$V9,";"),function(x) x[4])

# Create separate column for ensembl gene symbols
Ensembl$gene.sym<-ens.gene.name
Ensembl<-Ensembl[- grep("NT_", Ensembl$V1),]
Ensembl<-Ensembl[- grep("MT", Ensembl$V1),]
Ensembl$gene.sym<-sub("gene_name ","",Ensembl$gene.sym)

# Generate a GRanges object from gtf and nongtf data frames above for working with chr ranges
ensembl.ranges<-GRanges(seqnames=Rle(Ensembl$V1),ranges=IRanges(start=Ensembl$V4,end=Ensembl$V5),strand=Ensembl$V7,gene.names=Ensembl$gene.sym)
noncode.ranges<-GRanges(seqnames=Rle(Noncode$seqnames),ranges=IRanges(start=Noncode$start,end=Noncode$end),strand=Noncode$strand,gene.names=Noncode$gene.names)

# Combine both ensembl and noncode ranges (convert them to data frames first)
code.noncode <- rbind(as.data.frame(ensembl.ranges),as.data.frame(noncode.ranges))
code.noncode<-code.noncode[order(code.noncode$strand,code.noncode$seqnames,code.noncode$start),]

# Generate a GRanges object for the combined data frame
code.noncode.ranges<-GRanges(seqnames=Rle(code.noncode$seqnames),ranges=IRanges(start=code.noncode$start,end=code.noncode$end),strand=code.noncode$strand,gene.names=code.noncode$gene.names)


# Find overlaps between exons 
code.noncode.ovls<-findOverlaps(code.noncode.ranges,ignoreSelf=T,ignoreRedundant=T,type="any")

code.noncode.ovls.mat <- as.matrix(code.noncode.ovls)

code.noncode.ovls.dta <- data.frame(as.data.frame(code.noncode.ranges)[code.noncode.ovls.mat[,1],], as.data.frame(code.noncode.ranges)[code.noncode.ovls.mat[,2],])

# Verify if there are overlatps between exons of different genes
code.noncode.ovls.genes.inds <- which(code.noncode.ovls.dta$gene.names != code.noncode.ovls.dta$gene.names.1)

code.noncode.ovls.rm.inds <- as.integer(code.noncode.ovls.mat[code.noncode.ovls.genes.inds,])

# Prepare a list of the ranges in which exons are overlapping between genes
code.noncode.ovls.genes.df <- as.data.frame(code.noncode.ranges[unique(code.noncode.ovls.rm.inds)])

# Convert data frame to bed format so that we can run bedtools (intersectBed and subtractBed)
code.noncode.ovls.genes.df <-code.noncode.ovls.genes.df[c(1,2,3,6,4,5)]
code.noncode.ovls.genes.df$width<-0 #This should actually be the "score" column

# Convert whole list of exons to bed format (this is the list that bedtools will compare the above one against)
code.noncode<-code.noncode[c(1,2,3,6,4,5)]
code.noncode$width<-0

# Generate bed files for exons overlapping between genes and for the original list of exons 
write.table(code.noncode.ovls.genes.df,"CodeNoncode_overlapping_exons_btwn_genes.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(code.noncode,"CodeNoncode_exons.bed",quote=F,col.names=F,row.names=F,sep="\t")

# This was run to remove from original list just the portions of exons that overlap between genes
# sort -k1,1 -k2,2n overlapping_exons_btwn_genes.bed > overlapping_exons_btwn_genes_sorted.bed
# sort -k1,1 -k2,2n exons.bed > exons_sorted.bed
# /lawrencedata/darakjia/bedtools-2.17.0/bin/intersectBed -a exons_sorted.bed -b overlapping_exons_btwn_genes_sorted.bed > intersect.bed
# /lawrencedata/darakjia/bedtools-2.17.0/bin/subtractBed -s -a exons_sorted.bed -b intersect.bed > non_overlapping_exons.bed

# Get back the file with exons regions that do not overlap
no_gene_ovl_exons <- read.table("CodeNoncode_non_overlapping_exons.bed")

# Convert it to genomic ranges
no.gene.ovl.exons.ranges<-GRanges(seqnames=Rle(no_gene_ovl_exons$V1),ranges=IRanges(start=no_gene_ovl_exons$V2,end=no_gene_ovl_exons$V3),strand=no_gene_ovl_exons$V6,gene.names=no_gene_ovl_exons$V4)

# Merge exon isoforms
no.gene.ovl.exons.ranges.red<-reduce(split(no.gene.ovl.exons.ranges,elementMetadata(no.gene.ovl.exons.ranges)$gene.names))

# Convert above objects into data frame
no.gene.ovl.exons.ranges.red.df<-as.data.frame(no.gene.ovl.exons.ranges.red)

# Re-arrange data frame in bed format and write object to file
no.gene.ovl.exons.ranges.red.df<-no.gene.ovl.exons.ranges.red.df[c(3,4,5,2,6,7)]
no.gene.ovl.exons.ranges.red.df$width<-0

write.table(no.gene.ovl.exons.ranges.red.df,"CodeNoncode_reduced_non_overlapping_exons.bed", sep="\t",quote=F,col.names=F,row.names=F)

# RUN THIS OUTSIDE R
#################################################################################################
# Run coverage with bedtools
#################################################################################################
# Sort the loci bed file:
# sort -k6 -k1,1 -k2,2n CodeNoncode_reduced_non_overlapping_exons.bed > CodeNoncode_reduced_non_overlapping_exons_sorted.bed
# Add "chr" as a prefix to the chromosome name
# sed "s/^/chr/g" CodeNoncode_reduced_non_overlapping_exons_sorted.bed > CodeNoncode_reduced_non_overlapping_exons_sorted_chradd.bed
# Run the coverage script
# nohup ./run_bedtools_coverage_1.sh > run_bedtools_coverage_1.log 2>&1 </dev/null &
#################################################################################################

