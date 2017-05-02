################################################################################
##
## This script creates a table of statistics pertaining to the differential 
## expression of all mouse genes listed in the Mouse Genome Index. Here it
## gathers t-statistics, but this script was also used to create tables of
## log2(fold-changes) and multiple hypothesis testing corrected P-values. For 
## each microarry or RNA-seq experiment, it requires a table with columns for 
## MGI gene symbol and t-statistic (e.g. Limma's topTable output). Input 
## variables are capitalised. 
##
################################################################################
#
########################################
#
## Read in the results tables from the
## differential expression analysis of
## your microarray/RNA-seq experiments.
#
########################################
results <- c("results1","results2","results3,...,resultsN")
for(i in 1:length(results)){
tmp <- read.table(file=results[i],head=T,sep="\t",comment.char="",quote="")
assign(results[i],tmp)
}
########################################
#
## Create a list of the unique gene
## symbols included in the differential
## expression experiments that are also
## included in the Mouse Genome Index.
#
########################################
genes <- as.list(rep(NA,length(results)))
for(i in 1:length(results)){
j <- get(as.character(results[i]))
genes[[i]] <- as.character(j$SYMBOL)
}
genes <- as.character(unlist(genes))
genes.unique <- setdiff(genes,c(NA,""))
mgi.symbols <- read.table(file="PATH/TO/GENESYMBOLS.txt",head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F) 
genes.unique <- intersect(genes.unique,mgi.symbols$Symbol)
genes.unique <- sort(genes.unique)
########################################
#
## Create a table of t-statistcs for
## each gene in genes.unique.
#
########################################
t.table <- array(NA,dim=c(length(genes.unique),length(results)))
rownames(t.table) <- genes.unique
colnames(t.table) <- results
for(i in 1:length(results)){
j <- get(as.character(results[i]))
thisStudy <- data.frame(j$SYMBOL,(if("t" %in% colnames(j)) j$t else j$adj.P.Val))
colnames(thisStudy) <- c("SYMBOL","t")
mapped <- which(unique(thisStudy$SYMBOL) %in% genes.unique)
tofill <- which(genes.unique %in% thisStudy$SYMBOL)
new.statistics <- thisStudy$t[!duplicated(thisStudy$SYMBOL)][mapped]
t.table[tofill,i] <- new.statistics[order(unique(thisStudy$SYMBOL)[mapped])]
}
########################################
#
## Quantile normalise the table of t-
## statistics.
#
########################################
library(preprocessCore)
t.table.norm <- normalize.quantiles(t.table)
colnames(t.table.norm) <- colnames(t.table)
rownames(t.table.norm) <- rownames(t.table)
boxplot(t.table.norm)
