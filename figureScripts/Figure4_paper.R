# Michelle Webb
# Figure 4
# Spatial Transcriptomics ESTIMATE Analysis
# June 15, 2023
# Paths have been changed

library(Seurat) # version 4.3.0
library(ggplot2) # version 3.4.2
library(dplyr) # version 1.1.2
library(tLOH) # version 1.6.0

# Variables
distinctColors <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2",
                    "#FF00FF","palegreen2","darkturquoise","brown","#CAB2D6","#FDBF6F")

# Figure 4A
seurat <- Load10X_Spatial(data.dir = '~/DiffuseMidlineGlioma7/outs')
cluster <- read.csv('~/DiffuseMidlineGlioma7/outs/analysis/clustering/graphclust/clusters.csv', 
                    row.names = 1)
seurat@meta.data <- transform(merge(seurat@meta.data,cluster,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
SpatialDimPlot(seurat, group.by='Cluster', image.alpha=0,stroke=0, pt.size.factor = 0) +
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold')) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    labs(fill='Cluster')

# Figure 4B
clusterPlot <- SpatialDimPlot(seurat, group.by='Cluster', image.alpha=0,stroke=0) + 
    scale_fill_manual(values=distinctColors[1:9]) + 
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold')) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    labs(fill='Cluster')
ggsave(clusterPlot,filename = '~/diffuseMidlineGlioma7_clusterPlot.png')


# Figure 4C
estimateOutput <- read.csv('~/DiffuseMidlineGlioma7_estimateOutput.csv', row.names = 1)
seurat@meta.data <- transform(merge(seurat@meta.data,estimateOutput[c('StromalScore','TumorPurity')],by=0,all=TRUE),row.names=Row.names,Row.names=NULL)

tumorPurityPlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'TumorPurity', image.alpha=0) + theme(
    legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.position='right')
ggsave(tumorPurityPlot,
       filename = '~/diffuseMidlineGlioma7_tumorPurityPlot.png')

# Figure 4D
stromalScorePlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'StromalScore', image.alpha=0) + theme(
    legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.position='right')
ggsave(stromalScorePlot,
       filename = '~/diffuseMidlineGlioma7_tumorPurityPlot.png')

# Figure 4E
pdgfraPlot <- SpatialFeaturePlot(seurat, features = 'PDGFRA', pt.size.factor = 0) + theme(
    legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.position='none')

ggsave(pdgfraPlot,filename = '~/diffuseMidlineGlioma7_PDGFRA.png')

# Figure 4F
minidf <- seurat@meta.data
colors <- distinctColors[1:length(unique(minidf$Cluster))]
boxplot <- ggplot(minidf,aes(x=as.factor(Cluster),y=StromalScore,color=as.factor(Cluster))) + 
    geom_boxplot() + 
    scale_color_manual(values = colors) +
    theme(axis.text.y = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.position='none') +
    xlab('Cluster') +
    ylab('Stromal Score')
ggsave(boxplot,filename = '~/diffuseMidlineGlioma7_stromalScore_boxplot.png')

# Figure 4G
tlohPlot <- ggplot(toPlot,aes(x=position,
                        y=alleleFraction,
                        color=segment_State)) + 
    geom_point(size = 0.25) +
    facet_grid(cluster~chromosome, switch = 'x', scales = 'free') +
    scale_color_manual(values = c('blue','goldenrod'), labels = c('Heterozygous','LOH')) +
    theme(strip.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
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
          aspect.ratio = 1) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    labs(color='State') +
    xlab('Chromosome') +
    ylab('Allele Fraction')
ggsave(tlohPlot, filename = '~/diffuseMidlineGlioma7_paper_tLOH.png')
