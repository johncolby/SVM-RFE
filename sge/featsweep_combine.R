#! /usr/bin/env Rscript

# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

source('../msvmRFE.R')

files = list.files(pattern='feat_')
nfeat = sort(as.numeric(gsub('feat_([^\\.]+).Rdata', '\\1', files)))

# Load/combine results from sub-jobs
featsweep = list()
for(i in nfeat){
    load(sprintf('feat_%i.Rdata', i))
    featsweep[i] = list(results.tmp)
}

# Save out the combined results
save(featsweep, file='featsweep.Rdata')

# Calculate the "chance" rate
load('input.Rdata')
no.info = 1 - max(prop.table(table(input[,1])))

# Generate figures
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info)
dev.off()
dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info, ylim=c(0,0.5))
dev.off()

# Clean up
system('rm jobarray* feat_*')
