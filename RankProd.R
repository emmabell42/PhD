################################################################################
##
## This script performs Rank Products on a table of differential expression 
## log2(fold-changes).
##
################################################################################
library(RankProd)
########################################
#
## Create an origin vector to identify 
## different experiments, and a class
## vector to identify different data
## classes.
#
########################################
origin <- colnames(logFC.table)
logFCs.cl <- rep(1,length(results))
########################################
#
## Run RPadvance. This will generate two
## tables where for each gene:
##	- Table1 gives the statistical
##  significance of downregulation;
##	- Table2 gives the statistical
##	significance of upregulation.
#
########################################
RP.adv.out <- RPadvance(logFC.table,logFCs.cl,origin,num.perm=100,logged=TRUE,gene.names=rownames(logFC.table),rand=123)
RP.genes <- topGene(RP.adv.out,num.gene=nrow(logFC.table),method="pfp",logged=TRUE,logbase=2,gene.names=rownames(logFC.table))
down.genes <- RP.genes$Table1[RP.genes$Table1$pfp <= 0.05,]
up.genes <- RP.genes$Table2[RP.genes$Table2$pfp <= 0.05,]
