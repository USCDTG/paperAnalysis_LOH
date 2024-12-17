library(dplyr) # version 1.1.4  
library(purrr) #  version 1.0.2

GlioblastomaA1 <- read.csv('~/GlioblastomaA1/outs/metrics_summary.csv',check.names=FALSE)
GlioblastomaA2 <- read.csv('~/GlioblastomaA2/outs/metrics_summary.csv',check.names=FALSE)
GlioblastomaA3 <- read.csv('~/GlioblastomaA3/outs/metrics_summary.csv',check.names=FALSE)
GlioblastomaA4 <- read.csv('~/GlioblastomaA4/outs/metrics_summary.csv',check.names=FALSE)
GlioblastomaA5 <- read.csv('~/GlioblastomaA5/outs/metrics_summary.csv',check.names=FALSE)
OligodendrogliomaA1 <- read.csv('~/OligodendrogliomaA1/outs/metrics_summary.csv',check.names=FALSE)
OligodendrogliomaA2 <- read.csv('~/OligodendrogliomaA2/outs/metrics_summary.csv',check.names=FALSE)
OligodendrogliomaA3 <- read.csv('~/OligodendrogliomaA3/outs/metrics_summary.csv',check.names=FALSE)
AstrocytomaA1 <- read.csv('~/AstrocytomaA1/outs/metrics_summary.csv',check.names=FALSE)
AstrocytomaA2 <- read.csv('~/AstrocytomaA2/outs/metrics_summary.csv',check.names=FALSE)
DiffuseMidlineGlomaA1 <- read.csv('~/DiffuseMidlineGliomaA1/outs/metrics_summary.csv',check.names=FALSE)

samples <- list(GlioblastomaA1,GlioblastomaA2,GlioblastomaA3,GlioblastomaA4,GlioblastomaA5,
                OligodendrogliomaA1,OligodendrogliomaA2,OligodendrogliomaA3,AstrocytomaA1,AstrocytomaA2,DiffuseMidlineGlomaA1)

dataset <- reduce(samples,full_join)
dataset$`Sample ID` <- c('GlioblastomaA1','GlioblastomaA2','GlioblastomaA3','GlioblastomaA4','GlioblastomaA5',
                         'OligodendrogliomaA1','OligodendrogliomaA2','OligodendrodendrogliomaA3','AstrocytomaA1','AstrocytomaA2','DiffuseMidlineGlomaA1')

finalTable <- dataset[c('Sample ID',
                        'Number of Spots Under Tissue',
                        'Number of Reads',
                        'Total Genes Detected',
                        'Median Genes per Spot',
                        'Sequencing Saturation',
                        'Median UMI Counts per Spot')]

names(finalTable) <- c('Sample','Spots Under Tissue', 'Total Reads (Mils)','Total Genes Detected','Median Genes Per Spot',
                       'Sequencing Saturation','Median UMI counts per spot')

write.csv(finalTable, file = 'Supplementary_Table_2.csv',row.names=FALSE)


