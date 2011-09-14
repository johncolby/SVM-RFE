# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

# Start up an R session in the SVM-RFE directory. Then work through these commands.

# Set up R environment
set.seed(12345)
library(e1071)
source('msvmRFE.R')
load('demo/input.Rdata')

# Take a look at the expected input structure
dim(input)
input[1:5,1:5]

# Basic usage
svmRFE(input, k=10, halve.above=100)

# Set up cross validation
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds
folds = lapply(1:nfold, function(x) which(folds == x))
folds

# Perform feature ranking on all training sets
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
length(results)
results

# Obtain top features across ALL folds
top.features = WriteFeatures(results, input, save=F)
head(top.features)

# Estimate generalization error using a varying number of top features
featsweep = lapply(1:5, FeatSweep.wrap, results, input)
featsweep

# Make plot
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info)
dev.off()