# Figure 4
# Paths have been changed

library(Seurat) #version 5.1.0   
library(ggplot2) #version 3.5.1    
library(dplyr) #version 1.1.4
library(tLOH) #version 1.6.0
library(viridis) #version 0.6.5

setwd('~/Figures/')
seurat <- Load10X_Spatial('~/DiffuseMidlineGlioma7/outs/')
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- NormalizeData(seurat, 
                    normalization.method = "LogNormalize", 
                    scale.factor = 10000)

# Variables
distinctColors <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2",
                    "#FF00FF","palegreen2","darkturquoise","brown","#CAB2D6","#FDBF6F")

# Figure 4A
estimateOutput <- read.csv('~/FFSlideED1_estimateOutput.csv', row.names = 1)
seurat@meta.data <- transform(merge(seurat@meta.data,estimateOutput[c('StromalScore','TumorPurity')],by=0,all=TRUE),row.names=Row.names,Row.names=NULL)
cluster <- read.csv('~/DiffuseMidlineGlioma7/outs/analysis/clustering/graphclust/clusters.csv',row.names=1) #10X graph clusters
seurat@meta.data <- transform(merge(seurat@meta.data,cluster,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
clusterPlot <- SpatialDimPlot(seurat, group.by='Cluster', image.alpha=0,stroke=0,pt.size.factor = 2.6) + 
    scale_fill_manual(values=distinctColors[1:9]) + 
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold')) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    labs(fill='Cluster')
ggsave(clusterPlot,filename = 'Figure4_a_clusterPlot.png')

# Figure 4B
minidf <- seurat@meta.data[c('Cluster','StromalScore')]

get_box_stats <- function(StromalScore, upper_limit = max(minidf$StromalScore) * 1.15) {
    return(data.frame(
        y = 0.95 * upper_limit,
        label = paste(
            "Count =", length(StromalScore), "\n",
            "Mean =", round(mean(StromalScore), 2), "\n",
            "Median =", round(median(StromalScore), 2), "\n"
        )
    ))
}

colors <- distinctColors[1:length(unique(minidf$Cluster))]
boxplot <- ggplot(minidf,aes(x=as.factor(Cluster),y=StromalScore,color=as.factor(Cluster))) + 
    stat_boxplot(geom='errorbar', linetype=1, width=0.5)+  #whiskers
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", size=2) +
    stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.5,size=2.6,color='black') +
    #stat_summary(fun.data = mean_se, geom = "errorbar",linetype='dashed',color='darkgrey',linewidth=0.3) +
    scale_color_manual(values = colors) +
    theme(axis.text.y = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.line = element_line(color='black'),
          axis.title.x = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
          legend.position='none') +
    scale_y_continuous(limits = c(-600, 2000), breaks = c(-500,0,500,1000,1500)) +
    xlab('Cluster') +
    ylab('Stromal Score')
    
ggsave(boxplot,filename = 'Figure4_b_stromalScore_boxplot.png',height = 7,width=7)

# Figure 4C
tumorPurityPlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'TumorPurity', image.alpha=0,pt.size.factor = 2.5) + 
    theme(
    legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.position='right')
ggsave(tumorPurityPlot,
       filename = 'Figure4_c_tumorPurity_spatialPlot.png',height=7,width=7)

# Figure 4D
stromalScorePlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'StromalScore', image.alpha=0,pt.size.factor = 2.5) + theme(
    legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
    legend.position='right')
ggsave(stromalScorePlot,
       filename = 'Figure4_d_stromalScore_spatialPlot.png',height=7,width=7)

# Figure 4E
plot1 <- SpatialFeaturePlot(seurat,
                            features = 'PDGFRA',
                            images='slice1',
                            pt.size.factor = 2.6,
                            slot = 'data',
                            image.alpha=0) +
    scale_fill_viridis(option='viridis') +
    theme(legend.position = 'right',
          legend.text= element_text(face='bold',family='Arial Narrow',size=15),
          legend.title= element_text(face='bold.italic',family='Arial Narrow',size=15))
ggsave(plot1, file = 'Fig4_e_PDGFRA.png',height=7,width=7)


# Figure 4F
# Generated using Loupe Browser

# Figure 4G
# Generated with per-cluster SNP allele count files
toPlot <- read.csv('FFSlideED1_tLOH_perCluster.csv')
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
ggsave(tlohPlot, filename = 'Figure4_g_diffuseMidlineGlioma7_tLOH.png')



