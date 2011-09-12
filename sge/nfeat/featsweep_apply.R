#! /usr/bin/env Rscript

# Input arguments:
# i

# Parse command line arguments    
args=(commandArgs(TRUE))   
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

# Enter the code you want to parallelize...
library('e1071')

source('../../msvmRFE.R')

load('../results.Rdata')
load('../input.Rdata')

results.tmp = FeatSweep.wrap(i, results, input)

# Save out the results from this sub-job
save(results.tmp, file=sprintf('feat_%i.Rdata', i))