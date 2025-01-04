# Supplementary Figure 7
# Paths have been changed
library(Seurat) #version 5.1.0   
library(ggplot2) #version 3.5.1    
library(dplyr) #version 1.1.4
library(viridis) #version 0.6.5

setwd('~/Figures/ESTIMATE')
estimateList <- list.files('~/ESTIMATE_sourceData',full.names = TRUE)
samples <- gsub('_.*','',basename(estimateList))


for(i in 1:length(samples)){
    estimateData <- read.csv(estimateList[[i]],row.names = 1)
    seurat <- Load10X_Spatial(data.dir = paste0('~/spatialData/',samples[[i]],'/outs'))
    seurat@meta.data <- transform(merge(seurat@meta.data,estimateData,by=0,all=TRUE),row.names=Row.names,Row.names=NULL)
    ## Stromal
    stromalPlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'StromalScore', image.alpha=0,pt.size.factor = 2.5) + 
        theme(
            legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.position='right') +
        scale_fill_viridis(option='turbo')
    
    ggsave(stromalPlot,
           filename = paste0('SupplementaryFigure7_',samples[[i]],'_stromalScore_spatialPlot.png'),
           height=7,width=7)
    
    ## Immune
    immunePlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'ImmuneScore', image.alpha=0,pt.size.factor = 2.5) + 
        theme(
            legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.position='right') +
        scale_fill_viridis(option='turbo')
    
    ggsave(immunePlot,
           filename = paste0('SupplementaryFigure7_',samples[[i]],'_immuneScore_spatialPlot.png'),
           height=7,width=7)
    
    ## TumorPurity
    tumorPurityPlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'TumorPurity', image.alpha=0,pt.size.factor = 2.5) + 
        theme(
            legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.position='right') +
        scale_fill_viridis(option='turbo')
    
    ggsave(tumorPurityPlot,
           filename = paste0('SupplementaryFigure7_',samples[[i]],'_tumorPurity_spatialPlot.png'),
           height=7,width=7)
    
    ## ESTIMATE
    estimatePlot <- SpatialFeaturePlot(seurat, stroke=0,features = 'ESTIMATEScore', image.alpha=0,pt.size.factor = 2.5) + 
        theme(
            legend.text = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.title = element_text(size = 20, family = 'Arial Narrow', face = 'bold'),
            legend.position='right') +
        scale_fill_viridis(option='turbo')
    
    ggsave(estimatePlot,
           filename = paste0('SupplementaryFigure7_',samples[[i]],'_ESTIMATEScore_spatialPlot.png'),
           height=7,width=7)
}
