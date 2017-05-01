################################################################################
##
## This script loads Illumina microarray data into R, log2 transforms, quantile
## normalises, performs linear regression, and calls differentially expressed
## genes between two conditions. Input variables are capitalised.
##
################################################################################
#
########################################
#
## Load the relevant packages into R
#
########################################
#
library(lumi)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)
library(limma)
########################################
#
## Load the data into R for processing
## as a lumi batch object
#
########################################
data <- read.table(file="PATH/TO/FILE",head=T)
lumi <- lumiR(file="PATH/TO/FILE",lib.mapping="lumiMouseIDMapping")
lumi.exprs <- exprs(lumi)
lumi.exprs[which(lumi.exprs<1)] <- 1
########################################
#
## Log2 transform and quantile
## normalise the probe intensity values
#
########################################
log <- log2(lumi.exprs)
norm <- lumiN(log,method="quantile")
########################################
#
## Perform clustering analysis to
## check the arrays cluster as expected
## and boxplot to confirm normalised
## distribution
#
########################################
plot(hclust(dist(t(norm))))
boxplot(log(normexprs,base=2)
########################################
#
## Fit the probe intensities to a linear
## model
#
########################################
design <- array(0,dim=c(X,Y))
colnames(design) <- c("COLNAMES")
design[X,Y] <- 1
cont <- makeContrasts(CONTRAST="CONDITION1-CONDITION2",levels=design)
fit <- lmFit(norm,design)
fit <- contrasts.fit(fit,cont)
fit <- eBayes(fit)
########################################
#
## Annotate the probes with the names of
## their target genes
#
########################################
probe.list <- rownames(exprs(lumi))
if (require(lumiMouseAll.db) & require(annotate)) {
gene.symbol <- getSYMBOL(probe.list, 'lumiMouseAll.db')
fit$genes <- data.frame(ID=probe.list, SYMBOL=gene.symbol)
}
fit$genes$SYMBOL <- pData(featureData(data))$SYMBOL
########################################
#
## Output a table statistically
## evaluating differential expression
## between CONDITION1 and CONDITION2 
#
########################################
hits <- topTable(fit,number=nrow(norm))
