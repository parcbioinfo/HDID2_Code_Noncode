########################################################################################
# This script removes noncode rows present in the NONCODE2016_mouse_mm10_lncRNA.gtf file
# that have Ensembl counterparts listed in the Mus_musculus.GRCm38.78.gtf file under
# an actual gene name
#
# 09/13/2015
# Priscila Darakjian - OHSU - BNE
########################################################################################

library(GenomicRanges)
Ensembl.gtf<-read.table("gtfs/Mus_musculus.NCBIM37.67_exons.gtf",sep="\t", header=F,stringsAsFactors=F)
Noncode.gtf<-read.table("gtfs/NONCODEv4_mouse_lncRNA_exons.gtf",sep="\t", header=F,stringsAsFactors=F)

# Filter Ensembl.gtf to only other than protein-coding genes since
# those are the ones we are looking for redundancy with the Noncodes
Ensembl.gtf<-Ensembl.gtf[Ensembl.gtf$V2 %in% c("3prime_overlapping_ncrna","lincRNA","miRNA","ncrna_host","non_coding","snoRNA","snRNA","transcribed_unprocessed_pseudogene","unprocessed_pseudogene"),]
sym.gene.name.code<-sapply(strsplit(Ensembl.gtf$V9,";"),function(x) x[4])
sym.gene.name.noncode<-sapply(strsplit(Noncode.gtf$V9,";"),function(x) x[1])

Ensembl.gtf$gene.sym<-sym.gene.name.code
Noncode.gtf$gene.sym<-sym.gene.name.noncode

# For the Ensemble file Remove chromosomes named NT_... and MT, and 
# remove "gene_name" string from gene.names column
Ensembl.gtf<-Ensembl.gtf[- grep("NT_", Ensembl.gtf$V1),]
Ensembl.gtf$gene.sym<-sub("gene_name ","",Ensembl.gtf$gene.sym)

# For the NONCODE file remove "gene_name" string from gene.names column, and remove rows 
# with non-standard chroomosome names; also remove "chr" from the chromosome name since 
# Ensembl doesn't have that.
Noncode.gtf$gene.sym<-sub("gene_id ","",Noncode.gtf$gene.sym)
Noncode.gtf<-Noncode.gtf[-grep("chr13_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chr17_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chr1_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chr4_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chr8_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chr9_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chrM", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chrUn_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chrX_random", Noncode.gtf$V1),]
Noncode.gtf<-Noncode.gtf[-grep("chrY_random", Noncode.gtf$V1),]
Noncode.gtf$V1<-sub("chr","",Noncode.gtf$V1)

# Create Genomic Ranges from the gtf files
code.test.ranges<-GRanges(seqnames=Rle(Ensembl.gtf$V1),ranges=IRanges(start=Ensembl.gtf$V4,end=Ensembl.gtf$V5),strand=Ensembl.gtf$V7,gene.names=Ensembl.gtf$gene.sym)

noncode.test.ranges<-GRanges(seqnames=Rle(Noncode.gtf$V1),ranges=IRanges(start=Noncode.gtf$V4,end=Noncode.gtf$V5),strand=Noncode.gtf$V7,gene.names=Noncode.gtf$gene.sym)

noncode.test.ranges<-unique(noncode.test.ranges)
code.test.ranges<-unique(code.test.ranges)

#####################################################################################
# Find Exact Overlaps - Code (just a precaution) (there should be none)
#####################################################################################

code.ranges.ovls<-findOverlaps(code.test.ranges,ignoreSelf=T,ignoreRedundant=T,type="equal")
code.ranges.ovls.mat <- unique(as.matrix(code.ranges.ovls))

code.ranges.ovls.dta <- unique(data.frame(as.data.frame(code.test.ranges)[code.ranges.ovls.mat[,1],], as.data.frame(code.test.ranges)[code.ranges.ovls.mat[,2],]))


#####################################################################################
# Find Exact Overlaps - Noncode (just a precaution) (there should be none)
#####################################################################################

noncode.ranges.ovls<-findOverlaps(noncode.test.ranges,ignoreSelf=T,ignoreRedundant=T,type="equal")
noncode.ranges.ovls.mat <- unique(as.matrix(noncode.ranges.ovls))

noncode.ranges.ovls.dta <- unique(data.frame(as.data.frame(noncode.test.ranges)[noncode.ranges.ovls.mat[,1],], as.data.frame(noncode.test.ranges)[noncode.ranges.ovls.mat[,2],]))

#####################################################################################
# If the code and noncode overlaps above comes up empty, continue below
# i.e. find exact overlaps between the code and noncode ranges (we want to eliminate
# the noncodes that actually have gene symbols in Ensembl)
#####################################################################################

code.noncode<-rbind(as.data.frame(code.test.ranges),as.data.frame(noncode.test.ranges))
code.noncode<-code.noncode[order(code.noncode$strand,code.noncode$seqnames,code.noncode$start),]

code.noncode.ranges<-GRanges(seqnames=Rle(code.noncode$seqnames),ranges=IRanges(start=code.noncode$start,end=code.noncode$end),strand=code.noncode$strand,gene.names=code.noncode$gene.names)

code.noncode.equal.ovls<- findOverlaps(code.noncode.ranges,ignoreSelf=T,ignoreRedundant=T,type="equal")
code.noncode.equal.ovls.mat<-as.matrix(code.noncode.equal.ovls)
code.noncode.equal.ovls.dta<-data.frame(as.data.frame(code.noncode.ranges)[code.noncode.equal.ovls.mat[,1],],as.data.frame(code.noncode.ranges)[code.noncode.equal.ovls.mat[,2],])

# Verify that there are only NONCODEs in gene.names.1
code.noncode.equal.ovls.dta[- grep("NONMMUG",code.noncode.equal.ovls.dta$gene.names.1),] # should return no rows

# Fetch the indexes of subject hits (noncodes) so we can remove those redundant noncodes from the data frame
code.noncode.equal.ovls.rm <- as.integer(code.noncode.equal.ovls.mat[,2])

# Remove redundant noncodes from the combined data frame
code.noncode.cln <- as.data.frame(code.noncode.ranges)[- code.noncode.equal.ovls.rm,]

# Extract only the noncodes from the combined data frame and write the resulting df to a file
noncode.cln<-code.noncode.cln[grep("NONMMUG",code.noncode.cln$gene.names),]
write.table(noncode.cln,"NONCODEv4_mouse_lncRNA_exons_nonredundant.txt",sep="\t", row.names=F,col.names=T,quote=F)

save(noncode.cln,file="noncode.cln.RData")
