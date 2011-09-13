#! /usr/bin/env Rscript

# Input arguments:
# nfold

# Example:
# msvmRFE_submit.R nfold=10

# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

# Parse command line arguments    
args=(commandArgs(TRUE))   
for(i in 1:length(args)){
	eval(parse(text=args[i]))
}

# Determine number of rows
load('input.Rdata')
nrows = nrow(input)

# Save out fold ids
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
save(folds, file='folds.Rdata')

# Create SGE job array file
fid = file('jobarray.txt', 'w')
for(ifold in 1:nfold){
	write('msvmRFE_apply.R', fid, append=T)
}
close(fid)

# Submit job array to the SGE queue
jobID = system('fsl_sub -T 200 -l . -t jobarray.txt', intern=T)
system(paste('echo', jobID))

# Submit job to combine smaller sub-jobs
system(sprintf('fsl_sub -T 5 -l . -j %s msvmRFE_combine.R', jobID))

# Print job info to stdout
system(paste('echo Number of folds:', nfold))
