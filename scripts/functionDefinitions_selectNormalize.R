prepareCountsData = function (rawCountsFile,sampleNamesFile)
{
  #read raw data (counts generated with bedtools multiBamCov)
  geneReadsRaw=read.table(rawCountsFile,header=F,sep="\t",stringsAsFactors = F)
  
  # Change column 4 name to "gene_sym"
  names(geneReadsRaw)[4]<-"gene_sym"
  
  # Combine the chromosome number, start location, and gene symbol to create a unique id  
  # column for each exon
  geneReadsRaw$exon<-paste(geneReadsRaw$V1,geneReadsRaw$V2,geneReadsRaw$V3,geneReadsRaw$gene_sym,sep="_")
  
  # Create a data frame with gene symbol and exon read counts
  exon_counts<-geneReadsRaw[,7:ncol(geneReadsRaw)]
  exon_counts<-cbind(geneReadsRaw$gene_sym,exon_counts)
  
  # Calculate the total counts for each gene for each sample
  gene_counts<-ddply(exon_counts, 1, numcolwise(sum))
  
  # Change the row names of the exon data frame to the exon unique ids (created above)
  rownames(exon_counts)<-exon_counts$exon
  
  # Remove the gene symbol and exon id columns from the exon data frame
  exon_counts<-exon_counts[,2:ncol(exon_counts)-1]
  
  # Change the row names of the gene data frame to the gene symbols
  names(gene_counts)[1]<-"gene_sym"
  rownames(gene_counts)<-gene_counts$gene_sym
  
  # Remove the gene symbol column from the gene data frame
  gene_counts<-gene_counts[,2:ncol(gene_counts)]
  
  # Load sample names in the order they are in the coverage file
  samples <- read.table(sampleNamesFile,header=F,stringsAsFactors = F)
  
  names(gene_counts)<-samples[,]
  names(exon_counts)<-samples[,]
  
  save(samples, gene_counts, exon_counts,file="data/gene_and_exon_counts.RData")
  
}
 

normalizeCounts = function (exon_counts,gene_counts) 
{
  UQnormFactors_exons=calcNormFactors(exon_counts, method=c("upperquartile"))
  UQnormFactors_genes=calcNormFactors(gene_counts, method=c("upperquartile"))
  
  effectiveLibrarySizes_exons= UQnormFactors_exons*colSums(exon_counts)
  effectiveLibrarySizes_genes= UQnormFactors_genes*colSums(gene_counts)
  
  meanEffLibSize_exons=mean(effectiveLibrarySizes_exons)
  meanEffLibSize_genes=mean(effectiveLibrarySizes_genes)
  
  countNormFactor_exons= meanEffLibSize_exons/effectiveLibrarySizes_exons
  countNormFactor_genes= meanEffLibSize_genes/effectiveLibrarySizes_genes
  
  normalizedGeneCountsUQ_exons=0* exon_counts
  normalizedGeneCountsUQ_genes=0* gene_counts
  
  for (sample in 1:dim(normalizedGeneCountsUQ_exons)[2]){  
    normalizedGeneCountsUQ_exons[,sample]= exon_counts[, sample]* countNormFactor_exons[sample]	
  }
  
  for (sample in 1:dim(normalizedGeneCountsUQ_genes)[2]){  
    normalizedGeneCountsUQ_genes[,sample]= gene_counts[, sample]* countNormFactor_genes[sample]  
  }
  
  save(normalizedGeneCountsUQ_exons,normalizedGeneCountsUQ_genes,file="data/normalizedCounts.RData")
}

detectOutliers = function (gene_counts)
{
  geneMatrixData=as.matrix(gene_counts)
  colnames(geneMatrixData)=colnames(gene_counts)
  for (r in 1:dim(geneMatrixData)[1])
  {
    for (c in 1:dim(geneMatrixData)[2])
    {
      geneMatrixData[r,c]=gene_counts[r,c]
    }
  }
  geneIAC=cor(geneMatrixData,method="p") 
  
  # IAC shoud ideally be 0.8 or more - lower IACs could be outliers
  geneCluster=hclust(as.dist(1-geneIAC),method="average") 
  
  geneMeanIAC=apply(geneIAC,2,mean) 
  geneSdCorr=sd(geneMeanIAC) 
  geneNumberSd=(geneMeanIAC-mean(geneMeanIAC))/geneSdCorr 

  # samples more than 2 standard deviations from mean Inter Array Correlation (IAC) are considered outliers
  sdout=-2 
  geneOutliers=dimnames(geneMatrixData)[[2]][geneNumberSd<sdout] 
  save(geneMatrixData,geneIAC,geneCluster,geneMeanIAC,geneSdCorr,geneNumberSd, geneOutliers,file="data/outlierSamples.RData")
}

diffExpressionTwoGroups = function (gene_counts_group1,gene_counts_group2,groupNames)
{
  gn1<-groupNames[1]
  gn2<-groupNames[2]
  groupSelection=c(rep(gn1,dim(gene_counts_group1)[2]),rep(gn2,dim(gene_counts_group2)[2]))
  groupSelection=factor(groupSelection)
  
  d=DGEList(counts= cbind(gene_counts_group1, gene_counts_group2), group= groupSelection)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  de.tgw <- exactTest(d, dispersion="tagwise") 
  de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust.method="fdr")
  
  return(list(d,de.tgw,de.calls))
  
}

runDEAnalysis = function (de)
{
  # Multiple comparison correction
  # results from sgof come out sorted but un-named (!!!!) so the original pvalues and geneNames need to be sorted
  pValDE=de[[2]]$table$PValue
  geneNamesDE<-rownames(de[[2]]$table)
  names(pValDE)=geneNamesDE
  
  sortIndexes=sort.int(pValDE, decreasing = F, index.return=T)$ix
  sortedGeneNamesDE=geneNamesDE[sortIndexes]
  
  adjustedResultsDE<-SGoF(u=pValDE)
  summary(adjustedResultsDE)
  
  sortedAdjustedPvalsDE=adjustedResultsDE$Adjusted.pvalues
  names(sortedAdjustedPvalsDE)=sortedGeneNamesDE
  
  return(sortedAdjustedPvalsDE)
}

runDVAnalysis = function (geneNamesDE, normCountsG1,normCountsG2)
{
  pvalVar=rep(1, length(geneNamesDE))

  names(pvalVar)=geneNamesDE

  for (gene in geneNamesDE) 
  {
    pvalVar[gene]=var.test(x=t(normCountsG1[gene,]), y=t(normCountsG2[gene,]))$p.value
  }
  
  pvalVar[is.na(pvalVar)]=1

  sortIndexes=sort.int(pvalVar, decreasing = F, index.return=T)$ix
  sortedGeneNames=geneNamesDE[sortIndexes]
  adjustedResultsVar<-SGoF(u=pvalVar)
  sortedAdjustedPvalsDV=adjustedResultsVar$Adjusted.pvalues
  names(sortedAdjustedPvalsDV)=sortedGeneNames
  return(sortedAdjustedPvalsDV)
}