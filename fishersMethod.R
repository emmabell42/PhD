################################################################################
##
## This script generates a tabe of Fisher's combined P-values.
##
################################################################################
#
########################################
#
## Create the fishersMethod function
#
########################################
fishersMethod <- function(x) ifelse(sum(!is.na(x))>0,pchisq(-2 * sum(log(x[!is.na(x)])),df=2*length(x[!is.na(x)]),lower=FALSE),NA)
########################################
#
## Create a table to populate with
## the fishersMethod output and fill it.
#
########################################
fishers.table <- array(NA,dim=c(nrow(p.table.norm),1))
for(i in 1:nrow(p.table.norm)){
fishers.table[i,] <- fishersMethod(p.table.norm[i,])
}
rownames(fishers.table) <- rownames(p.table.norm)
colnames(fishers.table) <- "Combined P-value"
########################################
#
## FDR correct the Fisher's P-values and
## add these to the table.
#
########################################
table.p.adj <- as.matrix(p.adjust(fishers.table,method="fdr"))
fishers.table <- cbind(fishers.table,table.p.adj)
colnames(fishers.table) <- c("Combined P-value","Adj P-value")
fishers.table.order <- fishers.table[order(fishers.table[,2]),]
