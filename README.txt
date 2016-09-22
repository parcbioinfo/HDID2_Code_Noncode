


Here is the command used to get the count of genes for each biotype in the gtf file:
(in R)
typeCounts <- aggregate(Ensembl.gtf$gene_name, by=list("biotype"=Ensembl.gtf$gene_biotype), function(x) length(unique(x)))

