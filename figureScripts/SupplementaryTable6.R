# Supplementary Table 6
# Paths have been changed
library(dplyr) # version 1.1.4  
library(purrr) #  version 1.0.2

B6 <- read.csv('~/FFPE/SPA1_D/outs/metrics_summary.csv', check.names=FALSE)
B7 <- read.csv('~/FFPE/SPA2_A/outs/metrics_summary.csv', check.names=FALSE)
B8 <- read.csv('~/FFPE/SPA3_D/outs/metrics_summary.csv', check.names=FALSE)
B9 <- read.csv('~/FFPE/SPA5_D/outs/metrics_summary.csv', check.names=FALSE)
B10 <- read.csv('~/FFPE/SPA7_D/outs/metrics_summary.csv', check.names=FALSE)
B11 <- read.csv('~/FFPE/SPA8_A/outs/metrics_summary.csv', check.names=FALSE)

samples <- list(B6,B7,B8,B9,B10,B11)
dataset <- reduce(samples,full_join)
dataset$`Sample ID` <- c('Glioblastoma B6',
                       'Glioblastoma B7',
                       'Glioblastoma B8',
                       'Glioblastoma B9',
                       'Glioblastoma B10',
                       'Glioblastoma B11')

finalTable <- dataset[c('Sample ID',
          'Number of Spots Under Tissue',
          'Number of Reads',
          'Mean Reads per Spot',
          'Mean Reads Under Tissue per Spot',
          'Fraction of Spots Under Tissue',
          'Valid Barcodes',
          'Sequencing Saturation',
          'Median Genes per Spot',
          'Median UMI Counts per Spot')]

write.csv(finalTable, file = '~/Supplementary_Table_5.csv',row.names=FALSE)
