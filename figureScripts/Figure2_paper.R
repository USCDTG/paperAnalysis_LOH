# Michelle Webb
# Figure 2 Analysis Workflow
# June 15, 2023
# Paths have been changed

library(ggpubr) # version 0.6.0
library(GenomicRanges) # version 1.50.2
library(IRanges) # version 2.32.0
library(VariantAnnotation) # version 1.44.1
library(regioneR) # 1.30.0

# Figure 2A
# Diagram

# Figure 2B 
# IGV Image

# Figure 2C
options(scipen=999)
dataToUse <- read.csv('~/allSamples_tLOH.csv')
imported2 <- dataToUse %>% group_by(chromosome) %>%
    distinct(chromosome,position)
snps <- imported2[!(imported2$chromosome == 6
                    & imported2$position > 28510120 
                    & imported2$position < 33500500),]
snps$seqnames <- paste0('chr',snps$chromosome)
snps$chromosome <- NULL
snps$start <- snps$position
snps$end <- snps$position
snps$position <- NULL
snps$CHR <- gsub('^','chr',snps$chromosome)
snps <- snps %>% filter(seqnames !='chrX' & seqnames !='chrY')
rangedSNPS <- regioneR::toGRanges(A = data.frame(snps))


kp <- karyoploteR::plotKaryotype(genome="hg19")
kp <- karyoploteR::kpPlotDensity(kp, data=rangedSNPS, col="blue", window.size = 10000)
max.density <- kp$latest.plot$computed.values$max.density
# karyoploteR::kpAxis(kp, ymin=0, ymax=max.density, numticks = 3,cex =0.5)

png(filename = '~/Karyoplote.png', width = 520, height = 520)
pp <- karyoploteR::getDefaultPlotParams(plot.type=1)
pp$margin <- 0.1
pp$bottommargin <- 150
pp$data1height <- 500
pp$ideogramheight <- 50
kp <- karyoploteR::plotKaryotype(plot.type=1, genome='hg19',cex=1.85, plot.params = pp, chromosomes = c('chr1','chr2','chr3','chr4','chr5',
                                                                                                        'chr6','chr7','chr8','chr9','chr10',
                                                                                                        'chr11','chr12','chr13','chr14','chr15',
                                                                                                        'chr17','chr18','chr19','chr20','chr21',
                                                                                                        'chr22'))
kp <- karyoploteR::kpPlotDensity(kp, data=rangedSNPS, col = 'blue')
dev.off()
# Figure 2D
# Diagram

# Figure 2E
a1 <- ggplot(dataToUse,aes(x=segment_lengthOfInterval)) + geom_histogram(color = 'darkslategray4', fill = 'darkslategray4') +
    scale_x_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6),expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('Segment Length') +
    theme(axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'))

b1 <- ggplot(dataToUse,aes(x=segment_nSNPs)) + geom_histogram(color = 'darkslategray4', fill = 'darkslategray4') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('# of SNPs per Segment') +
    theme(axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'))

c1 <- ggplot(dataToUse,aes(x=segment_modePeak)) + geom_histogram(color = 'darkslategray4', fill = 'darkslategray4') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('Segment Mode Peak') +
    theme(axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'))
d1 <- ggplot(dataToUse,aes(x=segment_sequentialSum)) + geom_histogram(color = 'darkslategray4', fill = 'darkslategray4') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab('Segment Sequential Sum Values') +
    theme(axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'))

m <- ggpubr::ggarrange(a1,b1,c1,d1,
                       labels = c("i","ii",'iii','iv'),
                       ncol=2,
                       nrow=2)
ggsave(m,filename = '~/figure2_histograms.png')
