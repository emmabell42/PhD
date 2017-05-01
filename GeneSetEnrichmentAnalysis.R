################################################################################
##
## This script performs Gene Set Enrichment Analysis of candidate genes in 
## ConsensusPathDB pathways. It requires a list of genes with an associated
## statistic and the ConsensusPathDB list of biological pathways with their
## genes identified by gene symbol
##
################################################################################
#
########################################
#
## Read in ConsensusPathDB pathways and
## split the column containing the genes
## into a list
#
########################################
pathways <- read.table(file="PATH/TO/CPDB_pathways_genes.tab",sep="\t",header=TRUE)
pathway.genes <- list()
for(i in 1:nrow(pathways)){
        pathway.genes[[i]] <- strsplit(as.character(pathways$hgnc_gene_ids[i]),split=",")[[1]]
}
########################################
#
## Create a results table, populate it
## with the geneSetTest output and
## correct the P-values for multiple
## hypothesis testing
#
########################################
pathway.enrichments <- array(NA,dim=c(nrow(pathways),1))
rownames(pathway.enrichments) <- pathways$pathway
library(limma)
for(j in 1:nrow(pathways)){
		# default is for t-statistics. if using f-statistics (as in Fisher's meta-analysis statistics), comment out line 20 and uncomment line 21
                pathway.enrichments[j,1] <- geneSetTest(index=which(GENES[,3] %in% pathway.genes[[j]]),statistics=STATISTICS,ranks.only=FALSE,type="t",alternative="up")
		#pathway.enrichments[j,i] <- geneSetTest(selected=which(rownames(statistics.table) %in% pathway.genes[[j]]),statistics=statistics.table[,i],type="f",ranks.only=FALSE)
}
pathway.enrichments[,1] <- p.adjust(pathway.enrichments[,1],method="BH")
