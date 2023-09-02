# Michelle Webb
# Figure 3
# tLOH Analysis
# June 15, 2023

library(tLOH) # version 1.6.0
library(dplyr) # version 1.1.2
library(ggplot2) # version 3.4.2

# Variables
distinctColors <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2",
                    "#FF00FF","palegreen2","darkturquoise","brown","#CAB2D6","#FDBF6F")

# Function
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

# Figure 3A
seurat <- Load10X_Spatial(data.dir = '~/RecurrentGlioblastoma1/outs/')
cluster <- read.csv('~/RecurrentGlioblastoma1/analysis/clustering/graphclust/clusters.csv', row.names = 1)
seurat@meta.data <- transform(merge(seurat@meta.data,cluster,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

# Figure 3Ai
a <- SpatialDimPlot(seurat, group.by='cluster', image.alpha=0,stroke=0) + 
    scale_fill_manual(values=distinctColors[1:9]) + 
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold')) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    labs(fill='Cluster')
ggsave(a,filename = '~/RecurrentGlioblastoma1_clusterMap.png')


# Figure 3Aii
a <- SpatialFeaturePlot(seurat, features = 'EGFR', image.alpha = 0) + 
    theme(legend.position = 'right',
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'))
ggsave(a,filename = '~/RecurrentGlioblastoma1_EGFR.png')


# Figure 3Bi
seurat <- Load10X_Spatial(data.dir = '~/Glioblastoma10/outs/')

cluster <- read.csv('~/Glioblastoma10/outs/analysis/clustering/kmeans_3_clusters/clusters.csv', row.names = 1)

seurat@meta.data <- transform(merge(seurat@meta.data,cluster,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
a <- SpatialDimPlot(seurat, group.by='Cluster', image.alpha=0,stroke=0) + scale_fill_manual(values=c25[1:3]) + theme(legend.key = element_rect(colour = NA, fill = NA),
                                                                                                                     legend.box.background = element_blank(),
                                                                                                                     legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
                                                                                                                     legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold')) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    labs(fill='Cluster')
ggsave(a,filename = '~/Glioblastoma10_cluster.png')

# Figure 3Bii
a <- SpatialFeaturePlot(seurat, features = 'EGFR', image.alpha = 0) + theme(legend.position = 'right',
                                                                            legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
                                                                            legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'))
ggsave(a,filename = '~/Glioblastoma10_EGFR.png')


# Figure 3Ci

# Figure 3Cii

# Figure 3Ciii
a <- read.csv('~/KGLIUSCDC_0001_ps202012161514/cna/KGLIUSCDC_0001_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00013-KGLIUSCDC_0001_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00001_exo/KGLIUSCDC_0001_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00013-KGLIUSCDC_0001_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00001_exo.baf.txt', sep='\t')
b <- ggplot() +
    geom_rect(data=data.frame(Chromosome = c(1,3,5,9,11,13,15,17), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray82", alpha=0.5) +
    geom_rect(data=data.frame(Chromosome = c(7,10,19,21), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="mistyrose", alpha=0.5) +
    geom_point(size = 0.5,data = a, mapping = aes(x=Position,y=BAF,color = ifelse((Chromosome == '9' & Position < 43000000) | Chromosome == '7' | Chromosome == '10' | Chromosome == '19' | Chromosome == '20','goldenrod','blue'))) +
    # geom_point(data = a, mapping = aes(x=Position,y=BAF), size = 0.75) + 
    scale_color_manual(values = c('blue','goldenrod'), labels = c('Heterozygous','LOH')) +
    facet_grid(~Chromosome, switch = 'x',scales='free') +
    xlab('Chromosome') +
    ylab('B Allele Frequency') +
    theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.position = 'none',
          legend.background = element_rect(color = NA),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          panel.spacing = unit(0.00000001, 'lines'),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          aspect.ratio=5)
ggsave(b,filename = '~/RecurrentGlioblastoma1_exomeBAF.png')


# Figure 3Civ
a <- read.csv('~/KGLIUSCDC_0001_ps202012161514/cna/KGLIUSCDC_0001_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00013-KGLIUSCDC_0001_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00001_exo/KGLIUSCDC_0001_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00013-KGLIUSCDC_0001_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00001_exo.cna.tsv', sep='\t')
a <- a %>% filter(Chr != '24' & Chr != '23')
b <- ggplot() +
    geom_rect(data=data.frame(Chr = c(1,3,5,9,11,13,15,17), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray82", alpha=0.5) +
    geom_rect(data=data.frame(Chr = c(7,10,19,21), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="mistyrose", alpha=0.5) +
    geom_point(size = 0.25,data = a, mapping = aes(x=Position,y=Fold.Change,color = ifelse((Chr == '9' & Position < 43000000) | Chr == '7' | Chr == '10' | Chr == '19' | Chr == '20','goldenrod','blue'))) +
    # geom_point(data = a, mapping = aes(x=Position,y=BAF), size = 0.75) + 
    scale_color_manual(values = c('blue','goldenrod'), labels = c('Heterozygous','LOH')) +
    facet_grid(~Chr, switch = 'x',scales='free') +
    xlab('Chromosome') +
    ylab('Log2 Fold Change') +
    theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.position = 'none',
          legend.background = element_rect(color = NA),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          panel.spacing = unit(0.00000001, 'lines'),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          aspect.ratio=5)
ggsave(b,filename = '~/RecurrentGlioblastoma1_exomeCNV.png')


# Figure 3Di
sample <- read.csv('~/Glioblastoma10_alleleCounts.csv',check.names=FALSE)
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

# Figure 3Dii
ac <- list.files('~/Glioblastoma10/',full.names=TRUE)
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
ggsave(a,filename = '~/Glioblastoma10_tLOH_figure3.png')


data <- read.csv('~/Glioblastoma10_tLOH_perCluster.csv')

a <- plotSegments(data, c(7,9,10,13,15,19))

#Figure 3Diii
a <- read.csv('~/KGLIUSCDC_0011_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00023-KGLIUSCDC_0011_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00011_exo.baf.txt', sep='\t')
a <- a %>% filter(Chromosome != '24' & Chromosome != '23')
b <- ggplot() +
    geom_rect(data=data.frame(Chromosome = c(1,3,5,9,11,17,21), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray82", alpha=0.5) +
    geom_rect(data=data.frame(Chromosome = c(7,10,13,15,19), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="mistyrose", alpha=0.5) +
    geom_point(size = 0.5,data = a, mapping = aes(x=Position,y=BAF,color = ifelse((Chromosome == '9' & Position < 43000000) | Chromosome == '7' | Chromosome == '10' | (Chromosome == '13' & Position < 60000000) | (Chromosome == '15' & Position < 50000000) | Chromosome == '19','goldenrod','blue'))) +
    scale_color_manual(values = c('blue','goldenrod'), labels = c('Heterozygous','LOH')) +
    facet_grid(~Chromosome, switch = 'x',scales='free') +
    xlab('Chromosome') +
    ylab('B Allele Frequency') +
    theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.position = 'none',
          legend.background = element_rect(color = NA),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          panel.spacing = unit(0.00000001, 'lines'),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          aspect.ratio=5)
ggsave(b,filename = '~/Glioblastoma10_exomeBAF.png')


# Figure3Div
a <- read.csv('~/KGLIUSCDC_0011_01_WBLCONTTD_R00001S8A1M0001P0000_C1_AV6UOS_A00023-KGLIUSCDC_0011_01_GLITUM1TD_R00001S8A1M0001P0000_T1_AV6UOS_A00011_exo.cna.tsv', sep='\t')
a <- a %>% filter(Chr != '24' & Chr != '23')
b <- ggplot() +
    geom_rect(data=data.frame(Chr = c(1,3,5,9,11,17,21), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="gray82", alpha=0.5) +
    geom_rect(data=data.frame(Chr = c(7,10,13,15,19), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill="mistyrose", alpha=0.5) +
    geom_point(size = 0.25,data = a, mapping = aes(x=Position,y=Fold.Change,color = ifelse((Chr == '9' & Position < 43000000) | Chr == '7' | Chr == '10' | (Chr == '13' & Position < 60000000) | (Chr == '15' & Position < 50000000) | Chr == '19','goldenrod','blue'))) +
    # geom_point(data = a, mapping = aes(x=Position,y=BAF), size = 0.75) + 
    scale_color_manual(values = c('blue','goldenrod'), labels = c('Heterozygous','LOH')) +
    facet_grid(~Chr, switch = 'x',scales='free') +
    xlab('Chromosome') +
    ylab('Log2 Fold Change') +
    theme(strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.position = 'none',
          legend.background = element_rect(color = NA),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          panel.spacing = unit(0.00000001, 'lines'),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7, family = 'Arial Narrow', face = 'bold'),
          aspect.ratio=5)

ggsave(b,filename = '~/Glioblastoma10_exomeCNA.png')


