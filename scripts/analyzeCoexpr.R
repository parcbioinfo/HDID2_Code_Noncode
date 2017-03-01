library(foreach)
library(doMC)
registerDoMC()
library(multtest)
library(WGCNA)
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")
library(vegan)
library(ncf)

library(lawstat)

#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
#biocLite("GOstats")

getDoParWorkers()
options(cores=6)
getDoParWorkers()

#enableWGCNAThreads(nThreads = 27)

setwd("/lawrencedata/ongoing_analyses/RNASeq016/RNASeq016_noncode/NonCode_Analysis_2016")

source("scripts/functionDefinitions.R")

load("RNASeq016_code_noncode_sgof_adjusted_data.RData")
load("RNASeq016_code_noncode_DEDV_SignifData.RData")
geneNames_F=geneNames_DE_Fem_sel
geneNames_M=geneNames_DE_Mal_sel
save(geneNames_M,geneNames_F,file="geneNames.RData")
# divide the data in different groups
F_HDID<-HDIDFem_gene_counts_clean_norm[geneNames_F,]
M_HDID<-HDIDMal_gene_counts_clean_norm[geneNames_M,]
F_NPT<-NPTFem_gene_counts_clean_norm[geneNames_F,]
M_NPT<-NPTMal_gene_counts_clean_norm[geneNames_M,]
selectedGeneCounts_M<-cbind(M_HDID,M_NPT)
selectedGeneCounts_F<-cbind(F_HDID,F_NPT)

########################################################################################################################################
# For Females
adjConsensus_F=adjacency(t(selectedGeneCounts_F), power=1)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus_F, powerVector = powers, verbose = 5, moreNetworkConcepts=T)
plotNetConstruction(sft)
softPowerCoexpr_F=7
adjCoexpr_F=adjConsensus_F^softPowerCoexpr_F
adjCoexpr_F[is.na(adjCoexpr_F)]=0
connCoexpr_F=rowSums(adjCoexpr_F)
hist(connCoexpr_F, 100)
hierADJCoexpr_F = hclust(as.dist(1-adjCoexpr_F),method="average");
# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr_F=cutreeHybrid(dendro = hierADJCoexpr_F, distM=1-adjCoexpr_F, cutHeight = 0.9999, minClusterSize = 130, deepSplit = 4, verbose = 2, indent = 0)
colorsCoexpr_F = labels2colors(hybridCoexpr_F$labels)
names(colorsCoexpr_F)=geneNames_F
table(colorsCoexpr_F)
length(table(colorsCoexpr_F))
modulesCoexpr_F=names(table(colorsCoexpr_F))
sum(colorsCoexpr_F=="grey")

fileConnSummary<-file("data/resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Number of genes selected for network construction in females: ", length(geneNames_F), sep=''), fileConnSummary)
writeLines(paste("Number modules: ", length(table(colorsCoexpr_F)), sep=''), fileConnSummary)
writeLines(paste("Number grey genes:  ",  sum(colorsCoexpr_F=="grey"), sep=''), fileConnSummary)
close(fileConnSummary)

# For Males
adjConsensus_M=adjacency(t(selectedGeneCounts_M), power=1)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus_M, powerVector = powers, verbose = 5, moreNetworkConcepts=T)
plotNetConstruction(sft)
softPowerCoexpr_M=7
adjCoexpr_M=adjConsensus_M^softPowerCoexpr_M
adjCoexpr_M[is.na(adjCoexpr_M)]=0
connCoexpr_M=rowSums(adjCoexpr_M)
hist(connCoexpr_M, 100)
hierADJCoexpr_M = hclust(as.dist(1-adjCoexpr_M),method="average");
# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr_M=cutreeHybrid(dendro = hierADJCoexpr_M, distM=1-adjCoexpr_M, cutHeight = 0.999, minClusterSize = 130, deepSplit = 4, verbose = 2, indent = 0)
colorsCoexpr_M = labels2colors(hybridCoexpr_M$labels)
names(colorsCoexpr_M)=geneNames_M
table(colorsCoexpr_M)
length(table(colorsCoexpr_M))
modulesCoexpr_M=names(table(colorsCoexpr_M))
sum(colorsCoexpr_M=="grey")

fileConnSummary<-file("data/resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Number of genes selected for network construction in males: ", length(geneNames_M), sep=''), fileConnSummary)
writeLines(paste("Number modules: ", length(table(colorsCoexpr_M)), sep=''), fileConnSummary)
writeLines(paste("Number grey genes:  ",  sum(colorsCoexpr_M=="grey"), sep=''), fileConnSummary)
close(fileConnSummary)


adj_F_HDID=adjacency(t(F_HDID), power=softPowerCoexpr_F)
adj_F_NPT=adjacency(t(F_NPT), power=softPowerCoexpr_F)

adj_M_HDID=adjacency(t(M_HDID), power=softPowerCoexpr_M)
adj_M_NPT=adjacency(t(M_NPT), power=softPowerCoexpr_M)

save(F_HDID,F_NPT,M_NPT,M_HDID,selectedGeneCounts_F,selectedGeneCounts_M,adjConsensus_F,adjConsensus_M,adjCoexpr_F, adjCoexpr_M, hierADJCoexpr_M,hierADJCoexpr_F,hybridCoexpr_M,hybridCoexpr_F,colorsCoexpr_F,colorsCoexpr_M,modulesCoexpr_F, modulesCoexpr_M,adj_F_NPT,adj_F_HDID, adj_M_NPT,adj_M_HDID,file="data/adjCoexprModules.RData")


#load("data/adjCoexprModules.RData")
########################################################################################################
coexprConn_F=intramodularConnectivity(adjMat = adjCoexpr_F, colors = colorsCoexpr_F, scaleByMax = T)
coexprConnHDID_F=intramodularConnectivity(adjMat = adj_F_HDID, colors = colorsCoexpr_F, scaleByMax = F)
coexprConnNPT_F=intramodularConnectivity(adjMat = adj_F_NPT, colors = colorsCoexpr_F, scaleByMax = F)

coexprConn_M=intramodularConnectivity(adjMat = adjCoexpr_M, colors = colorsCoexpr_M, scaleByMax = T)
coexprConnHDID_M=intramodularConnectivity(adjMat = adj_M_HDID, colors = colorsCoexpr_M, scaleByMax = F)
coexprConnNPT_M=intramodularConnectivity(adjMat = adj_M_NPT, colors = colorsCoexpr_M, scaleByMax = F)

save(coexprConn_M,coexprConn_F,coexprConnNPT_M,coexprConnNPT_F,coexprConnHDID_F,coexprConnHDID_M, file="coexprConn.RData")
##############################################################
########################################################################################################
colorsCoexpr_F_gsym<-colorsCoexpr_F
colorsCoexpr_M_gsym<-colorsCoexpr_M

names(colorsCoexpr_F_gsym)<-sub(":.*", "", names(colorsCoexpr_F_gsym))
names(colorsCoexpr_M_gsym)<-sub(":.*", "", names(colorsCoexpr_M_gsym))

neuronsList=read.csv("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons_F = moduleEnrichment (colorsCoexpr_F_gsym, neuronsSymbols)
moduleEnrichmentAstros_F = moduleEnrichment (colorsCoexpr_F_gsym, astrosSymbols)
moduleEnrichmentOligos_F = moduleEnrichment (colorsCoexpr_F_gsym, oligosSymbols)

moduleEnrichmentNeurons_M = moduleEnrichment(colorsCoexpr_M_gsym, neuronsSymbols)
moduleEnrichmentAstros_M = moduleEnrichment (colorsCoexpr_M_gsym, astrosSymbols)
moduleEnrichmentOligos_M = moduleEnrichment (colorsCoexpr_M_gsym, oligosSymbols)

cellTypeEnrichment_F=round(cbind(moduleEnrichmentNeurons_F,moduleEnrichmentAstros_F, moduleEnrichmentOligos_F),4)
cellTypeEnrichment_M=round(cbind(moduleEnrichmentNeurons_M,moduleEnrichmentAstros_M, moduleEnrichmentOligos_M),4)

colnames(cellTypeEnrichment_F)=c("Neurons", "Astros", "Oligos")
colnames(cellTypeEnrichment_M)=c("Neurons", "Astros", "Oligos")

rownames(cellTypeEnrichment_F)=modulesCoexpr_F
rownames(cellTypeEnrichment_M)=modulesCoexpr_M

write.csv(cellTypeEnrichment_F, file="data/resultsCoexpr/cellTypeEnrich_F.csv")
write.csv(cellTypeEnrichment_M, file="data/resultsCoexpr/cellTypeEnrich_M.csv")

modulesEnrichedNeurons_F=modulesCoexpr_F[cellTypeEnrichment_F[,"Neurons"] < (0.05/length(modulesCoexpr_F))]
modulesEnrichedAstros_F=modulesCoexpr_F[cellTypeEnrichment_F[,"Astros"] < (0.05/length(modulesCoexpr_F))]
modulesEnrichedOligos_F=modulesCoexpr_F[cellTypeEnrichment_F[,"Oligos"] < (0.05/length(modulesCoexpr_F))]

modulesEnrichedNeurons_M=modulesCoexpr_M[cellTypeEnrichment_M[,"Neurons"] < (0.05/length(modulesCoexpr_M))]
modulesEnrichedAstros_M=modulesCoexpr_M[cellTypeEnrichment_M[,"Astros"] < (0.05/length(modulesCoexpr_M))]
modulesEnrichedOligos_M=modulesCoexpr_M[cellTypeEnrichment_M[,"Oligos"] < (0.05/length(modulesCoexpr_M))]

modulesNeuros_F=""
for(i in 1:length(modulesEnrichedNeurons_F)){
  modulesNeuros_F=paste(modulesNeuros_F, modulesEnrichedNeurons_F[i])
}

modulesNeuros_M=""
for(i in 1:length(modulesEnrichedNeurons_M)){
  modulesNeuros_M=paste(modulesNeuros_M, modulesEnrichedNeurons_M[i])
}

modulesAstros_F=""
for(i in 1:length(modulesEnrichedAstros_F)){
  modulesAstros_F=paste(modulesAstros_F, modulesEnrichedAstros_F[i])
}

modulesAstros_M=""
for(i in 1:length(modulesEnrichedAstros_M)){
  modulesAstros_M=paste(modulesAstros_M, modulesEnrichedAstros_M[i])
}

modulesOligos_F=""
for(i in 1:length(modulesEnrichedOligos_F)){
  modulesOligos_F=paste(modulesOligos_F, modulesEnrichedOligos_F[i])
}

modulesOligos_M=""
for(i in 1:length(modulesEnrichedOligos_M)){
  modulesOligos_M=paste(modulesOligos_M, modulesEnrichedOligos_M[i])
}

fileConnSummary<-file("data/resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Modules enriched in neuronal markers in females: ", modulesNeuros_F, sep=""), fileConnSummary)
writeLines(paste("Modules enriched in astrocyte markers in females: ", modulesAstros_F, sep=" "), fileConnSummary)
writeLines(paste("Modules enriched in oligodendrocyte markers in females: ", modulesOligos_F, sep=" "), fileConnSummary)
writeLines(paste("Modules enriched in neuronal markers in males: ", modulesNeuros_M, sep=""), fileConnSummary)
writeLines(paste("Modules enriched in astrocyte markers in males: ", modulesAstros_M, sep=" "), fileConnSummary)
writeLines(paste("Modules enriched in oligodendrocyte markers in males: ", modulesOligos_M, sep=" "), fileConnSummary)

close(fileConnSummary)



########################################################################################################

# save the results below for use with enrinchR
try(dir.create("data/resultsCoexpr/moduleGeneList"), silent = T)

coexprConnConsensus_F=intramodularConnectivity(adjCoexpr_F,  colorsCoexpr_F, scaleByMax=T)
coexprConnConsensus_M=intramodularConnectivity(adjCoexpr_M,  colorsCoexpr_M, scaleByMax=T)

#totalScaledConnectivity=coexprConnConsensus[,"kTotal"]/max(coexprConnConsensus[,"kTotal"])

coexprConnHDID_F=intramodularConnectivity(adj_F_HDID,  colorsCoexpr_F, scaleByMax=T)
coexprConnNPT_F=intramodularConnectivity(adj_F_NPT,  colorsCoexpr_F, scaleByMax=T)

coexprConnHDID_M=intramodularConnectivity(adj_M_HDID,  colorsCoexpr_M, scaleByMax=T)
coexprConnNPT_M=intramodularConnectivity(adj_M_NPT,  colorsCoexpr_M, scaleByMax=T)

save(coexprConn_M,coexprConn_F,coexprConnNPT_M,coexprConnHDID_F,coexprConnNPT_F,coexprConnHDID_M,coexprConnConsensus_M,coexprConnConsensus_F,file="coexprConn.RData")
coexprResultsTable_F=cbind(colorsCoexpr_F, round(coexprConnConsensus_F[,"kWithin"],3), round(coexprConnHDID_F[,"kWithin"],3), round(coexprConnNPT_F[,"kWithin"],3))
coexprResultsTable_M=cbind(colorsCoexpr_M, round(coexprConnConsensus_M[,"kWithin"],3), round(coexprConnHDID_M[,"kWithin"],3), round(coexprConnNPT_M[,"kWithin"],3))

colnames(coexprResultsTable_F)=c("module", "consensus conn", "HDID Conn", "NPT Conn")
rownames(coexprResultsTable_F)=geneNames_F

colnames(coexprResultsTable_M)=c("module", "consensus conn", "HDID Conn", "NPT Conn")
rownames(coexprResultsTable_M)=geneNames_M

for (module in modulesCoexpr_F){
  print(module)
  currModuleInfo=cbind(rownames(coexprConnConsensus_F)[colorsCoexpr_F==module],round(coexprConnConsensus_F[colorsCoexpr_F==module,"kWithin"],2))
  colnames(currModuleInfo)=c("gene name", "module connectivity")
  write.csv(currModuleInfo, file=paste("data/resultsCoexpr/moduleGeneList/module_F_", module, ".csv", sep=""), row.names=F, col.names=F)  
}

for (module in modulesCoexpr_M){
  print(module)
  currModuleInfo=cbind(rownames(coexprConnConsensus_M)[colorsCoexpr_M==module],round(coexprConnConsensus_M[colorsCoexpr_M==module,"kWithin"],2))
  colnames(currModuleInfo)=c("gene name", "module connectivity")
  write.csv(currModuleInfo, file=paste("data/resultsCoexpr/moduleGeneList/module_M_", module, ".csv", sep=""), row.names=F, col.names=F)  
}

#############################################################################
# GO annotations

load("/lawrencedata/ongoing_analyses/RNA160225TP/all_fastq/Alignment_mm10/Read_counts_and_normalization/data/transcriptInfoMouse.RData")
setwd("data/")

annotateMouseModulesGO(colorsCoexpr_F_gsym, transcriptInfoMouse, type="Coexpr")
file.rename("resultsCoexpr/modulesCoexprGO.csv", "resultsCoexpr/modulesCoexprGO_F.csv")

annotateMouseModulesGO(colorsCoexpr_M_gsym, transcriptInfoMouse, type="Coexpr")
file.rename("resultsCoexpr/modulesCoexprGO.csv", "resultsCoexpr/modulesCoexprGO_M.csv")


##############################################################################
# record differential expression
load("RNASeq016_code_noncode_DEDV_SignifData.RData")

intersect(geneNames_F, rownames(results_final_selected_genes_Fem))
setdiff(geneNames_F, rownames(results_final_selected_genes_Fem))
results_DEDVNet_F<-results_final_selected_genes_Fem[geneNames_F,]

intersect(geneNames_M, rownames(results_final_selected_genes_Mal))
setdiff(geneNames_M, rownames(results_final_selected_genes_Mal))
results_DEDVNet_M<-results_final_selected_genes_Mal[geneNames_M,]

save(results_DEDVNet_M,results_DEDVNet_F,file="Results_DEDVNet.RData")
##################################################################################33
# Evaluate changes in edge strength FOR FEMALES
##################################################################################
load("resultsCoexpr/adjCoexprModules.RData")
rawAdj_HDID_F=adjacency(t(F_HDID), power=1)
rawAdj_NPT_F=adjacency(t(F_NPT), power=1)

#######################################################################################
#The folowing two commands were run in EXACLOUD because of "lack of memory" in lawrence
#diffEdges_F = diffEdges(rawAdj_HDID_F, rawAdj_NPT_F, n1=dim(F_HDID)[2], n2=dim(F_NPT)[2], pThreshold=0.01, adjThreshold=0.5, nCores=7)
#save(diffEdges_F, file="resultsCoexpr/diffEdges_F.RData")

######################################
# Procedure for running in exacloud:
######################################
save(rawAdj_NPT_F,rawAdj_HDID_F,F_HDID,F_NPT,file="resultsCoexpr/rawAdj_counts_F.RData")
# Transfer required files to exacloud - 
# scp resultsCoexpr/rawAdj_counts_F.RData darakjia@exacloud.ohsu.edu:/home/exacloud/lustre1/PARC/darakjia/projects/interm
# scp ../scripts/functionalDefinitions.R darakjia@exacloud.ohsu.edu:/home/exacloud/lustre1/PARC/darakjia/projects/interm

# Create the following R script (run_diffEdges_F.R):
# source("functionDefinitions.R")
# library(foreach)
# library(doMC)
# registerDoMC()
# getDoParWorkers()
# options(cores=10)
# getDoParWorkers()
# load("rawAdj_counts_F.RData")
# 
# diffEdges_F = diffEdges(rawAdj_HDID_F, rawAdj_NPT_F, n1=dim(F_HDID)[2], n2=dim(F_NPT)[2], pThreshold=0.01, adjThreshold=0.5, nCores=10)
# 
# save(diffEdges_F, file="diffEdges_F.RData")

# Request a cluster for an interactive session:
# condor_submit -interactive -append 'request_cpus=10' -append 'request_memory = 120 GB'

# Run external R to process R script
# R CMD run_diffEdges_F.R 

load("diffEdges_F.RData")
load("geneNames.RData")

totalEdges_F=(length(geneNames_F))^2
affectedEdges_F=sum(diffEdges_F)
edgeChangeRate_F=affectedEdges_F/totalEdges_F

geneChangeEdgeCount_F=rowSums(diffEdges_F)
names(geneChangeEdgeCount_F)=geneNames_F

pValuesEdgeChange_F=rep(1,length(geneNames_F))
names(pValuesEdgeChange_F)=geneNames_F

for (gene in geneNames_F){
  pValuesEdgeChange_F[gene]=binom.test(x=geneChangeEdgeCount_F[gene], n=length(geneNames_F), p=edgeChangeRate_F, alternative  ="g")$p.value
}

pValues_F=pValuesEdgeChange_F
names(pValues_F)=geneNames_F

sortIndexes_F=sort.int(pValues_F, decreasing = F, index.return=T)$ix
sortedGeneNames_F=geneNames_F[sortIndexes_F]

adjustedResults_F<-SGoF(u=pValues_F)
summary(adjustedResults_F)

sortedAdjustedPvals_DW_F=adjustedResults_F$Adjusted.pvalues
names(sortedAdjustedPvals_DW_F)=sortedGeneNames_F


##################################################################################33
# Evaluate changes in edge strength FOR MALES
##################################################################################33

rawAdj_HDID_M=adjacency(t(M_HDID), power=1)
rawAdj_NPT_M=adjacency(t(M_NPT), power=1)

#######################################################################################
#The folowing two commands were run in EXACLOUD because of "lack of memory" in lawrence
#diffEdges_M = diffEdges(rawAdj_HDID_M, rawAdj_NPT_M, n1=dim(M_HDID)[2], n2=dim(M_NPT)[2], pThreshold=0.01, adjThreshold=0.5, nCores=7)
#save(diffEdges_M, file="resultsCoexpr/diffEdges_M.RData")

######################################
# Procedure for running in exacloud:
######################################
save(rawAdj_NPT_M,rawAdj_HDID_M,M_HDID,M_NPT,file="resultsCoexpr/rawAdj_counts_M.RData")
# Transfer required files to exacloud - 
# scp resultsCoexpr/rawAdj_counts_M.RData darakjia@exacloud.ohsu.edu:/home/exacloud/lustre1/PARC/darakjia/projects/interm
# scp ../scripts/functionalDefinitions.R darakjia@exacloud.ohsu.edu:/home/exacloud/lustre1/PARC/darakjia/projects/interm

# Create the following R script (run_diffEdges_M.R):
# source("functionDefinitions.R")
# library(foreach)
# library(doMC)
# registerDoMC()
# getDoParWorkers()
# options(cores=10)
# getDoParWorkers()
# load("rawAdj_counts_M.RData")
# 
# diffEdges_M = diffEdges(rawAdj_HDID_M, rawAdj_NPT_M, n1=dim(M_HDID)[2], n2=dim(M_NPT)[2], pThreshold=0.01, adjThreshold=0.5, nCores=10)
# 
# save(diffEdges_M, file="diffEdges_M.RData")

# Request a cluster for an interactive session:
# condor_submit -interactive -append 'request_cpus=10' -append 'request_memory = 120 GB'

# Run external R to process R script
# R CMD run_diffEdges_M.R 

load("diffEdges_M.RData")
load("geneNames.RData")

totalEdges_M=(length(geneNames_M))^2
affectedEdges_M=sum(diffEdges_M)
edgeChangeRate_M=affectedEdges_M/totalEdges_M

geneChangeEdgeCount_M=rowSums(diffEdges_M)
names(geneChangeEdgeCount_M)=geneNames_M

pValuesEdgeChange_M=rep(1,length(geneNames_M))
names(pValuesEdgeChange_M)=geneNames_M

for (gene in geneNames_M){
  pValuesEdgeChange_M[gene]=binom.test(x=geneChangeEdgeCount_M[gene], n=length(geneNames_M), p=edgeChangeRate_M, alternative  ="g")$p.value
}

pValues_M=pValuesEdgeChange_M
names(pValues_M)=geneNames_M

sortIndexes_M=sort.int(pValues_M, decreasing = F, index.return=T)$ix
sortedGeneNames_M=geneNames_M[sortIndexes_M]

adjustedResults_M<-SGoF(u=pValues_M)
summary(adjustedResults_F)

sortedAdjustedPvals_DW_M=adjustedResults_M$Adjusted.pvalues
names(sortedAdjustedPvals_DW_M)=sortedGeneNames_M

save(pValuesEdgeChange_M,pValuesEdgeChange_F,file="pValuesEdgeChange.RData")
save(sortedAdjustedPvals_DW_F,sortedAdjustedPvals_DW_M,file="sortedAdjustedPVals_DiffEdge.RData")
####################################################################################
# Summary Results for Females
####################################################################################
load("Results_DEDVNet.RData")
load("coexprConn.RData")
load("adjCoexprModules.RData")
load("diffEdges_F.RData")
load("sortedAdjustedPVals_DiffEdge.RData")
summaryResults_F=cbind(colorsCoexpr_F, results_DEDVNet_F[,c(1:2, 4:11)], round(coexprConnConsensus_F[,"kWithin"],3), round(coexprConnHDID_F[,"kWithin"],3), round(coexprConnNPT_F[,"kWithin"],3),round(pValuesEdgeChange_F,3), round(sortedAdjustedPvals_DW_F[geneNames_F],3),geneChangeEdgeCount_F)
colnames(summaryResults_F)
colnames(summaryResults_F)=c("module","logFC", "logCPM", "mean.counts.HDID","mean.counts.NPT","p.val.DE", "adj.p.DE" , "sd.counts.HDID","sd.counts.NPT","p.val.DV", "adj.p.DV","modular conn Consensus","modular conn HDID", "modular conn NPT","p.val.DiffEdge", "adj.p.DiffEdge", "changed edge count")

#sanity checks
plot(summaryResults_F[,"logFC"], -log10(summaryResults_F[,"p.val.DE"]))
plot(summaryResults_F[,"changed edge count"], -log10(summaryResults_F[,"p.val.DiffEdge"]))
plot(summaryResults_F[,"changed edge count"], -log10(summaryResults_F[,"adj.p.DiffEdge"]))

write.csv(summaryResults_F, file="resultsCoexpr/tableResults_F.csv")
save(summaryResults_F,file="summaryResults_F.RData")
####################################################################################
# Summary Results for Males
####################################################################################
load("Results_DEDVNet.RData")
load("coexprConn.RData")
load("adjCoexprModules.RData")
load("diffEdges_M.RData")
load("sortedAdjustedPVals_DiffEdge.RData")
summaryResults_M=cbind(colorsCoexpr_M, results_DEDVNet_M[,c(1:2, 4:11)], round(coexprConnConsensus_M[,"kWithin"],3), round(coexprConnHDID_M[,"kWithin"],3), round(coexprConnNPT_M[,"kWithin"],3),round(pValuesEdgeChange_M,3), round(sortedAdjustedPvals_DW_M[geneNames_M],3),geneChangeEdgeCount_M)
colnames(summaryResults_M)
colnames(summaryResults_M)=c("module","logFC", "logCPM", "mean.counts.HDID","mean.counts.NPT","p.val.DE", "adj.p.DE" , "sd.counts.HDID","sd.counts.NPT","p.val.DV", "adj.p.DV","modular conn Consensus","modular conn HDID", "modular conn NPT","p.val.DiffEdge", "adj.p.DiffEdge", "changed edge count")

#sanity checks
plot(summaryResults_M[,"logFC"], -log10(summaryResults_M[,"p.val.DE"]))
plot(summaryResults_M[,"changed edge count"], -log10(summaryResults_M[,"p.val.DiffEdge"]))
plot(summaryResults_M[,"changed edge count"], -log10(summaryResults_M[,"adj.p.DiffEdge"]))

write.csv(summaryResults_M, file="resultsCoexpr/tableResults_M.csv")
save(summaryResults_M,file="summaryResults_M.RData")
###################################################################################
# Summary Results for module enrichment for females
#################################################################################
deGenes_F=geneNames_F[summaryResults_F[,"p.val.DE"] < 0.01]
dvGenes_F=geneNames_F[summaryResults_F[,"p.val.DV"] < 0.01]
dwGenes_F=geneNames_F[summaryResults_F[,"p.val.DiffEdge"] < 0.01]

modulesEnrichDE_F=moduleEnrichment(colorsCoexpr_F, deGenes_F)
affectedModulesDE_F=names(modulesEnrichDE_F)[modulesEnrichDE_F<(0.05/length(modulesCoexpr_F))]

modulesDE_F=""
for (i in 1:length(affectedModulesDE_F)){
  modulesDE_F=paste(modulesDE_F, affectedModulesDE_F[i], sep=" ")
}

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Modules affected by DE expression changes in females: ", modulesDE_F, sep=' '), fileConnSummary)
close(fileConnSummary)

modulesEnrichDV_F=moduleEnrichment(colorsCoexpr_F, dvGenes_F)
affectedModulesDV_F=names(modulesEnrichDV_F)[modulesEnrichDV_F<(0.05/length(modulesCoexpr_F))]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

modulesDV_F=""
for (i in 1:length(affectedModulesDV_F)){
  modulesDV_F=paste(modulesDV_F, affectedModulesDV_F[i], sep=" ")
}

writeLines(paste("Modules affected by DV changes in females: ", modulesDV_F, sep=' '), fileConnSummary)
close(fileConnSummary)

modulesEnrichDW_F=moduleEnrichment(colorsCoexpr_F, dwGenes_F)
affectedModulesDW_F=names(modulesEnrichDW_F)[modulesEnrichDW_F<(0.05/length(modulesCoexpr_F))]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
modulesDW_F=""
for (i in 1:length(affectedModulesDW_F)){
  modulesDW_F=paste(modulesDW_F, affectedModulesDW_F[i], sep=" ")
}

writeLines(paste("Modules affected by DiffEdge changes in females: ", modulesDW_F, sep=' '), fileConnSummary)
close(fileConnSummary)

#####################################################################################
# Summary Results for module enrichment for males
#################################################################################
deGenes_M=geneNames_M[summaryResults_M[,"p.val.DE"] < 0.01]
dvGenes_M=geneNames_M[summaryResults_M[,"p.val.DV"] < 0.01]
dwGenes_M=geneNames_F[summaryResults_M[,"p.val.DiffEdge"] < 0.01]

modulesEnrichDE_M=moduleEnrichment(colorsCoexpr_M, deGenes_M)
affectedModulesDE_M=names(modulesEnrichDE_M)[modulesEnrichDE_M<(0.05/length(modulesCoexpr_M))]

modulesDE_M=""
for (i in 1:length(affectedModulesDE_M)){
  modulesDE_M=paste(modulesDE_M, affectedModulesDE_M[i], sep=" ")
}

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Modules affected by DE expression changes in males: ", modulesDE_M, sep=' '), fileConnSummary)
close(fileConnSummary)

modulesEnrichDV_M=moduleEnrichment(colorsCoexpr_M, dvGenes_M)
affectedModulesDV_M=names(modulesEnrichDV_M)[modulesEnrichDV_M<(0.05/length(modulesCoexpr_M))]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")

modulesDV_M=""
for (i in 1:length(affectedModulesDV_M)){
  modulesDV_M=paste(modulesDV_M, affectedModulesDV_M[i], sep=" ")
}

writeLines(paste("Modules affected by DV changes in males: ", modulesDV_M, sep=' '), fileConnSummary)
close(fileConnSummary)

modulesEnrichDW_M=moduleEnrichment(colorsCoexpr_M, dwGenes_M)
affectedModulesDW_M=names(modulesEnrichDW_M)[modulesEnrichDW_M<(0.05/length(modulesCoexpr_M))]

fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
modulesDW_M=""
for (i in 1:length(affectedModulesDW_M)){
  modulesDW_M=paste(modulesDW_M, affectedModulesDW_M[i], sep=" ")
}

writeLines(paste("Modules affected by DiffEdge changes in males: ", modulesDW_M, sep=' '), fileConnSummary)
close(fileConnSummary)
 
save(modulesEnrichDW_M,modulesEnrichDV_M,modulesEnrichDE_M,modulesEnrichDW_F,modulesEnrichDV_F,modulesEnrichDE_F,affectedModulesDW_M,affectedModulesDV_M,affectedModulesDE_M,affectedModulesDW_F,affectedModulesDV_F,affectedModulesDE_F, file="modulesEnrich.RData")
