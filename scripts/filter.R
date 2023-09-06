library(dplyr)
args = commandArgs(trailingOnly=TRUE)
path = args[1]

a <- list.files(path, full.names = TRUE, pattern='*csv')
b <- lapply(a,read.csv,check.names=FALSE)
c <- lapply(b, function(x) x %>% filter(grepl('UTR|missense|nonsynonymous',INFO)) %>%
                filter(!grepl('pseudogene|intron|intragenic',INFO)) %>%
                filter(QUAL > 500) %>%
                filter(grepl('protein_coding',INFO)) %>% 
                filter(grepl('rs',ID)))
for(i in 1:length(a)){
    a[[i]] <- gsub('.csv','_filtered.csv',a[[i]])
    write.csv(c[[i]], file = a[[i]], row.names=FALSE)
}
