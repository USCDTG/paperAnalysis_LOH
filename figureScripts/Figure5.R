library(tLOH) # version 1.6.0
library(dplyr) # version 1.1.2
library(ggplot2) # version 3.4.2

# Figure 5a.i-v
# GATK Somatic copy number variant discovery tool

# Figure 5b.i-v
# GATK Somatic copy number variant discovery tool

# Figure 5a.vi-vii
# Generated with Loupe Browser

# Figure 5b.vi-vii
# Generated with Loupe Browser

# Figure 5a.viii-x
plotSegments <- function(data, listToHighlight){
    default <- c(1,3,5,9,11,13,15,17,19,21)
    if(missing(listToHighlight)) {
        x <- default
        ggplot(data) + 
            geom_rect(data=data.frame(chromosome = x, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray86", alpha=0.5) +
            geom_segment(data=data,aes(x=segment_intervalStart,
                                       xend=segment_intervalEnd,
                                       y=0.5,yend=0.5,
                                       color=segment_State), size = 2) +
            facet_grid(cluster~chromosome, switch = 'x', scales = 'free') +
            scale_color_manual(values = c('blue','goldenrod','grey'), labels = c('Heterozygous','LOH','Undefined')) +
            xlab("Chromosome") +
            ylab("") +
            labs(color = c('State')) +
            theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.position = 'bottom',
                  legend.background = element_rect(color = NA),
                  legend.key = element_rect(colour = NA, fill = NA),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(0.00000001, 'lines'),
                  panel.grid = element_blank(),
                  panel.background = element_rect(fill = NA, color = NA),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  aspect.ratio = 0.75) +
            guides(colour = guide_legend(override.aes = list(size=3)))
    } else {
        x <- setdiff(default,listToHighlight)
        ggplot(data) + 
            geom_rect(data=data.frame(chromosome = x, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray86", alpha=0.5) +
            geom_rect(data=data.frame(chromosome = listToHighlight, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="mistyrose", alpha=0.5) +
            geom_segment(data=data,aes(x=segment_intervalStart,
                                       xend=segment_intervalEnd,
                                       y=0.5,yend=0.5,
                                       color=segment_State), size = 2) +
            facet_grid(cluster~chromosome, switch = 'x', scales = 'free') +
            scale_color_manual(values = c('blue','goldenrod','grey'), labels = c('Heterozygous','LOH','Undefined')) +
            xlab("Chromosome") +
            ylab("") +
            labs(color = c('State')) +
            theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
                  legend.position = 'bottom',
                  legend.background = element_rect(color = NA),
                  legend.key = element_rect(colour = NA, fill = NA),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(0.00000001, 'lines'),
                  panel.grid = element_blank(),
                  panel.background = element_rect(fill = NA, color = NA),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  aspect.ratio = 0.75) +
            guides(colour = guide_legend(override.aes = list(size=3)))
    }
}
importFromAC <- function(file){ # importing data from allele count files
    a <- read.csv(file,check.names = FALSE)
    b <- a[c('#CHROM','POS','rsID','REF_COUNT_CALC','ALT_COUNT_CALC','sample')]
    names(b) <- c('CHR','POS','rsID','REF','ALT','CLUSTER')
    b$TOTAL <- b$ALT + b$REF
    b$CHR <- gsub('chr','',b$CHR)
    return(b)
}
data("initialStartProbabilities")
trProbs <- cbind(c(0.999,0.001),c(0.001,0.999))
sample <- read.csv('~/GlioblastomaA1_alleleCounts.csv',check.names=FALSE)
sample <- sample %>% filter(QUAL > 500)
sample <- sample %>% filter(rsID != '.')
sample <- sample %>% filter(grepl('UTR|missense|nonsynonymous',INFO))

b <- sample[c('#CHROM','POS','rsID','REF_COUNT_CALC','ALT_COUNT_CALC','sample')]
names(b) <- c('CHR','POS','rsID','REF','ALT','sample')
b$CLUSTER <- 1
b$TOTAL <- b$ALT + b$REF
b$CHR <- gsub('chr','',b$CHR)


merged <- b
merged$CLUSTER <- gsub('cluster','',merged$CLUSTER)

calc1 <- tLOHCalcUpdate(merged, 1.25,1.25,500,500,4)
calc2 <- hiddenMarkovAnalysis(calc1, initialStartProbabilities, trProbs)
calc2$regionScoreThreshold <- (calc2$segment_nSNPs * 0.75) * 0.5
calc3 <- dplyr::mutate(.data = calc2, segment_State = dplyr::case_when(state == 'errored' ~ 'Errored',
                                                                       segment_sequentialSum > regionScoreThreshold ~ 'LOH',
                                                                       segment_sequentialSum > 20 & segment_modePeak > 0.01 & 
                                                                           (segment_medianBF < 0.5 & segment_medianBF > 0) ~ 'LOH',
                                                                       segment_modePeak > 0.01 & segment_medianBF > 1 ~ 'LOH',
                                                                       segment_medianBF < 0 ~ 'HET',
                                                                       TRUE ~ 'Undefined'))
a <- plotSegments(calc3, c(7,9,10,13,15,19))
ggsave(a,filename = 'Figure5_a_viii_GlioblastomaA1_tLOH.png')

ac <- list.files('~/GlioblastomaA1/',full.names=TRUE)
imported <- lapply(ac, function(x) importFromAC(x))
merged <- purrr::reduce(imported,dplyr::full_join)
merged$CLUSTER <- gsub('cluster','',merged$CLUSTER)
calc1 <- tLOHCalcUpdate(merged, 1.25,1.25,500,500,4)
calc2 <- hiddenMarkovAnalysis(calc1, initialStartProbabilities, trProbs)
calc2$regionScoreThreshold <- (calc2$segment_nSNPs * 0.75) * 0.5
calc3 <- dplyr::mutate(.data = calc2, segment_State = dplyr::case_when(state == 'errored' ~ 'Errored',
                                                                       segment_sequentialSum > regionScoreThreshold ~ 'LOH',
                                                                       segment_sequentialSum > 20 & segment_modePeak > 0.01 & 
                                                                           (segment_medianBF < 0.5 & segment_medianBF > 0) ~ 'LOH',
                                                                       segment_modePeak > 0.01 & segment_medianBF > 1 ~ 'LOH',
                                                                       segment_medianBF < 0 ~ 'HET',
                                                                       TRUE ~ 'Undefined'))
a <- plotSegments(calc3, c(7,10,13,15,19))
ggsave(a,filename = 'Figure5_a_ix_GlioblastomaA1_tLOH.png')

# Figure 5b.viii-x
sample <- read.csv('~/GlioblastomaA5_alleleCounts.csv',check.names=FALSE)
sample <- sample %>% filter(QUAL > 500)
sample <- sample %>% filter(rsID != '.')
sample <- sample %>% filter(grepl('UTR|missense|nonsynonymous',INFO))

b <- sample[c('#CHROM','POS','rsID','REF_COUNT_CALC','ALT_COUNT_CALC','sample')]
names(b) <- c('CHR','POS','rsID','REF','ALT','sample')
b$CLUSTER <- 1
b$TOTAL <- b$ALT + b$REF
b$CHR <- gsub('chr','',b$CHR)


merged <- b
merged$CLUSTER <- gsub('cluster','',merged$CLUSTER)

calc1 <- tLOHCalcUpdate(merged, 1.25,1.25,500,500,4)
calc2 <- hiddenMarkovAnalysis(calc1, initialStartProbabilities, trProbs)
calc2$regionScoreThreshold <- (calc2$segment_nSNPs * 0.75) * 0.5
calc3 <- dplyr::mutate(.data = calc2, segment_State = dplyr::case_when(state == 'errored' ~ 'Errored',
                                                                       segment_sequentialSum > regionScoreThreshold ~ 'LOH',
                                                                       segment_sequentialSum > 20 & segment_modePeak > 0.01 & 
                                                                           (segment_medianBF < 0.5 & segment_medianBF > 0) ~ 'LOH',
                                                                       segment_modePeak > 0.01 & segment_medianBF > 1 ~ 'LOH',
                                                                       segment_medianBF < 0 ~ 'HET',
                                                                       TRUE ~ 'Undefined'))
a <- plotSegments(calc3, c(7,9,10,13,15,19))
ggsave(a,filename = 'Figure5_b_viii_GlioblastomaA5_tLOH.png')

ac <- list.files('~/GlioblastomaA5/',full.names=TRUE)
imported <- lapply(ac, function(x) importFromAC(x))
merged <- purrr::reduce(imported,dplyr::full_join)
merged$CLUSTER <- gsub('cluster','',merged$CLUSTER)
calc1 <- tLOHCalcUpdate(merged, 1.25,1.25,500,500,4)
calc2 <- hiddenMarkovAnalysis(calc1, initialStartProbabilities, trProbs)
calc2$regionScoreThreshold <- (calc2$segment_nSNPs * 0.75) * 0.5
calc3 <- dplyr::mutate(.data = calc2, segment_State = dplyr::case_when(state == 'errored' ~ 'Errored',
                                                                       segment_sequentialSum > regionScoreThreshold ~ 'LOH',
                                                                       segment_sequentialSum > 20 & segment_modePeak > 0.01 & 
                                                                           (segment_medianBF < 0.5 & segment_medianBF > 0) ~ 'LOH',
                                                                       segment_modePeak > 0.01 & segment_medianBF > 1 ~ 'LOH',
                                                                       segment_medianBF < 0 ~ 'HET',
                                                                       TRUE ~ 'Undefined'))
a <- plotSegments(calc3, c(7,10,13,15,19))
ggsave(a,filename = 'Figure5_b_ix_GlioblastomaA5_tLOH.png')
