#! /usr/bin/env Rscript

ifold = as.numeric(Sys.getenv('SGE_TASK_ID'))

# Enter the code you want to parallelize...
library('e1071')

source('../msvmRFE.R')

load('folds.Rdata')
load('input.Rdata')

results.tmp = svmRFE.wrap(folds[[ifold]], input, k=10, halve.above=5000)

# Save out the results from this sub-job
save(results.tmp, file=sprintf('fold_%i.Rdata', ifold))
