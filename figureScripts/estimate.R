# ESTIMATE Analysis R Script
# .gct files were edited outside of R in excel
# paths renamed

library(estimate)

args = commandArgs(trailingOnly=TRUE)

sample = args[1]
inputFile = args[2]

filterCommonGenes(input.f = inputFile, output.f = paste0("~/estimateFINAL/",sample,'_filtered.gct'), id = "GeneSymbol")
estimateScore(input.ds = paste0("~/estimateFINAL/",sample,'_filtered.gct'), output.ds = paste0("~/estimateFINAL/",sample,'_estimate.gct'))


input <- read.csv(inputFile,sep='\t',check.names=FALSE,skip=2)
output <- read.csv(paste0("~/estimateFINAL/",sample,"_estimate.gct"),sep='\t',check.names=FALSE,skip=2)

input$GeneSymbol <- NULL
output[3] <- NULL
output[1] <- NULL
names(output) <- names(input)
write.csv(output,paste0("~/estimateFINAL/",sample,'_estimate_renamed.csv'),row.names=FALSE)

