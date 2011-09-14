#! /usr/bin/env Rscript

# Input arguments:
# nfeat

# Example:
# featsweep_submit.R nfeat=1:500
# featsweep_submit.R 'nfeat=c(1:100, seq(102, 500, by=2))'

# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

# Parse command line arguments    
args=(commandArgs(TRUE))   
for(i in 1:length(args)){
	eval(parse(text=args[i]))
}

# Create SGE job array file
fid = file('jobarray.txt', 'w')
for(i in nfeat){
    text = sprintf("featsweep_apply.R i=%i", i)
    write(text, fid, append=T)
}
close(fid)

# Submit job array to the SGE queue
jobID = system('fsl_sub -T 60 -l . -t jobarray.txt', intern=T)
system(sprintf('echo %s', jobID))

# Submit job to combine smaller sub-jobs
system(sprintf('fsl_sub -T 5 -l . -j %s featsweep_combine.R', jobID))

# Print job info to stdout
system(sprintf('echo Number of features to sweep: %i', length(nfeat)))
