library(Seurat)  #version 5.1.0 
library(dplyr) #version 1.1.4
library(ggplot2) #version 3.5.1

setwd('~/Figures')
# Data import
m <- readRDS('~/integrated_glioblastoma_samples_cohortB.rds')
m@meta.data$final_clusters <- as.factor(as.numeric(m@meta.data$SCT_snn_res.0.19))
m@meta.data <- m@meta.data %>% mutate(sample = recode(orig.ident,
                                                      'spa1' = 'Glioblastoma B6',
                                                      'spa2' = 'Glioblastoma B7',
                                                      'spa3' = 'Glioblastoma B8',
                                                      'spa5' = 'Glioblastoma B9',
                                                      'spa7' = 'Glioblastoma B10',
                                                      'spa8' = 'Glioblastoma B11'))
m@meta.data <- m@meta.data %>% mutate(cluster_annotation = recode(final_clusters,
                                                                  '1' = 'Stromal',
                                                                  '2' = 'OPC',
                                                                  '3' = 'Tumor \nAssociated \nMacrophages',
                                                                  '4' = 'Hypoxic-TAM',
                                                                  '5' = 'ECMs',
                                                                  '6' = 'Inflammation \nMicroglia',
                                                                  '7' = 'Hypoxic Invasive',
                                                                  '8' = 'Hypoxic MES',
                                                                  '9' = 'MES Astrocytes',
                                                                  '10' = 'Vascularization'))


# Figure 6a
# Plot was created with 

# Figure 6b
# Plot was created using Seurat DoHeatmap




# Figure 6c (1)
p <- DimPlot(m, reduction = 'umap',group.by='sample',pt.size = 0.5) +
    scale_color_manual(values = c('Glioblastoma B11' = '#E9CB58',
                                  'Glioblastoma B10' = '#67A055',
                                  'Glioblastoma B9' = '#82B5B2',
                                  'Glioblastoma B8' = '#D4605B',
                                  'Glioblastoma B7' = '#E78A3A',
                                  'Glioblastoma B6' = '#5577A5'),
                       limits = c("Glioblastoma B6", "Glioblastoma B7", "Glioblastoma B8",
                                  'Glioblastoma B9','Glioblastoma B10','Glioblastoma B11')) +
    ggtitle('') +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(0.1, 0.2),
          legend.text = element_text(face='bold',family='Arial Narrow',size=12)) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave(p, file = 'Fig6_c_sample_umap.png', height = 7, width = 7)

# Figure 6c (2)
DimPlot(m, reduction = 'umap',group.by='final_clusters',pt.size = 0.5, label=FALSE) +
    ggtitle('') +
    scale_color_manual(values = c('1' = '#9C9C9C',
                                  '2' = '#A6312B',
                                  '3' = '#851AF5',
                                  '4' = '#D0F163',
                                  '5' = '#20C1FF',
                                  '6' = '#EB873F',
                                  '7' = '#7CC8CE',
                                  '8' = '#377F80',
                                  '9' = '#81F490',
                                  '10' = '#3B3838')) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          # legend.position = c(0.1, 0.2),
          legend.position = 'none',
          legend.text = element_text(face='bold',family='Arial Narrow',size=12))
ggsave(p, file = 'Fig6_c_cluster_umap.png', height = 7, width = 7)

# Figure 6c (3)
stackedBarData <- m@meta.data %>% 
    group_by(sample,final_clusters) %>% 
    summarize(count = n())


p <- ggplot(stackedBarData, aes(fill=factor(sample,
                                       levels = c('Glioblastoma B11','Glioblastoma B10',
                                                  'Glioblastoma B9','Glioblastoma B8',
                                                  'Glioblastoma B7', 'Glioblastoma B6')), y=count, x=as.factor(final_clusters))) + 
    geom_bar(position="stack", stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color='black'),
          legend.text = element_text(size = 14,face='bold',family='Arial Narrow'),
          legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          strip.text = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.text.y = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          plot.title = element_text(size = 30, family = 'Arial Narrow', face = 'bold', hjust = 0.5)) +
    xlab('Cluster') +
    ylab('Number of Spots') +
    ggtitle('') +
    scale_fill_manual(values = c('Glioblastoma B11' = '#E9CB58',
                                  'Glioblastoma B10' = '#67A055',
                                  'Glioblastoma B9' = '#82B5B2',
                                  'Glioblastoma B8' = '#D4605B',
                                  'Glioblastoma B7' = '#E78A3A',
                                  'Glioblastoma B6' = '#5577A5'),
                      limits = c("Glioblastoma B6", "Glioblastoma B7", "Glioblastoma B8",
                                 'Glioblastoma B9','Glioblastoma B10','Glioblastoma B11')) +
    guides(size = guide_legend(override.aes = list(shape = 1)))
ggsave(p, file = 'Fig6_c_cluster_distribution_bar_plot.png', height = 7, width = 7)


# Figure d (1)
imageList <- c('slice1','slice1.2','slice1.3','slice1.4','slice1.5','slice1.6')
for(i in 1:length(imageList)){
    print(i)
    p <- SpatialDimPlot(m, 
                   group.by='final_clusters', 
                   images = imageList[[i]],
                   image.alpha=0,
                   pt.size.factor = 2.3) +
        scale_fill_manual(values = c('1' = '#9C9C9C',
                                     '2' = '#A6312B',
                                     '3' = '#851AF5',
                                     '4' = '#D0F163',
                                     '5' = '#20C1FF',
                                     '6' = '#EB873F',
                                     '7' = '#7CC8CE',
                                     '8' = '#377F80',
                                     '9' = '#81F490',
                                     '10' = '#3B3838')) +
        theme(legend.position = 'right') +
        ggtitle(unique(m@meta.data$sample)[[i]]) +
        guides(fill = guide_legend(override.aes = list(size = 5)))
    ggsave(p, file = paste0('Fig6_d_',unique(m@meta.data$sample)[[i]],'_spatialClusterMap.png'), height = 7, width = 7)
}

# Plots in e and f were generated using LoupeBrowser

