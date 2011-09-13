#! /usr/bin/env Rscript

# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

source('../msvmRFE.R')

# Load/combine results from sub-jobs
files = list.files(pattern='fold_')

results = list()
for(file in files){
    load(file)
    results = c(results, list(results.tmp))
}

# Save out a text file of the top features
load('input.Rdata')
WriteFeatures(results, input)

# Save out the combined results
save(results, file='results.Rdata')

# Clean up
system('rm jobarray* fold_* folds.Rdata')
