################################################################################
##
## This script loads Affymetrix microarray data into R, log2 transforms, 
## quantile normalises, performs linear regression, and calls differentially 
## expressed genes between two conditions. Input variables are capitalised.
##
################################################################################
#
########################################
#
## Load the relevant packages into R
#
########################################
#
library(affy)
library(limma)
########################################
#
## Load the data into R for processing
## as an Affy batch object
#
########################################
rawdata <- ReadAffy()
targets <- readTargets("TARGETS",row.names="filename")
exprs <- pm(rawdata)
########################################
#
## Log2 transform and quantile
## normalise the probe intensity values
#
########################################
normdata <- rma(rawdata)
normexprs <- exprs(normdata)
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
library(annotate)
library(ARRAY)
fit$genes$Symbol <- getSYMBOL(fit$genes$ID,"ARRAY")
########################################
#
## Output a table statistically
## evaluating differential expression
## between CONDITION1 and CONDITION2 
#
########################################
hits <- topTable(fit,number=nrow(norm))
