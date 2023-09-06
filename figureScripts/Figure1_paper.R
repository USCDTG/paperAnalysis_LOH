# Michelle Webb
# Glioma Dataset Integration
# June 15, 2023
# Paths have been changed

library(Seurat) # version 4.3.0
library(ggplot2) # version 3.4.2
library(tidyverse) # version 2.0.0
library(ComplexHeatmap) # version 2.14.0


# Variables
distinctColors <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2",
                    "#FF00FF","palegreen2","darkturquoise","brown","#CAB2D6","#FDBF6F")

# Figure 1A
samples <- c('FFSlideDB1', 'FFSlideDC1', 'FFSlideDD1', 'FFSlideEB1','FFSlideEC1',
             'FFSlideED1','FFSlideEE1','FFSlideGB1','FFSlideGD1','FFSlideGE1')
listOfSamples <- list()
for(i in 1:length(samples)){
    a <- Load10X_Spatial(paste0('~/spatial/',samples[[i]],'/outs'))
    a@meta.data$sample <- samples[[i]]
    a[["pctMito"]] <- PercentageFeatureSet(a, pattern = "^MT-")
    listOfSamples[[i]] <- a
}

additionalSample <- Load10X_Spatial(paste0('~/spatial/FF_D1'))
additionalSample@meta.data$sample <- 'FFD1'
additionalSample[["pctMito"]] <- PercentageFeatureSet(additionalSample, pattern = "^MT-")
listOfSamples2 <- unlist(list(additionalSample,listOfSamples))


intermediate <- lapply(listOfSamples2, SCTransform, assay='Spatial', method = "glmGamPoi", vars.to.regress = c('pctMito'))
features.list <- SelectIntegrationFeatures(object.list = intermediate)
integratedData <- lapply(intermediate, RunPCA, features=features.list)
integratedData <- PrepSCTIntegration(integratedData, anchor.features = features.list, assay = 'SCT')
anchors <- FindIntegrationAnchors(object.list = integratedData, anchor.features = features.list, 
                                          reduction = "rpca",   # reciprocal PCA
                                          normalization.method = 'SCT', 
                                          dims = 1:30)
integratedData <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "SCT")
integratedData <- RunPCA(integratedData, verbose = FALSE)
integratedData <- RunUMAP(integratedData, dims = 1:30)
# saveRDS(integratedData, '/Users/michelgw/integratedSampleDataset.rds')

sampleUMAP <- DimPlot(integratedData, reduction = 'umap', group.by = 'sample') +
    ggtitle('Integrated Samples') +
    scale_color_discrete(labels = c('Recurrent Glioblastoma 1 ','Oligodendroglioma 2 ','Glioblastoma 3','Recurrent Glioblastoma 4','Astrocytoma 5',
                                    'Astrocytoma 6','Diffuse Midline Glioma 7','Oligodendroglioma 8',
                                    'Glioblastoma 9','Glioblastoma 10','Oligodendroglioma 11')) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color='black'),
          strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold')) +
    guides(color = guide_legend(override.aes = list(size=5)))
ggsave(sampleUMAP, filename = 'integratedClusters_UMAP_highlightedSamples.png')

# Figure 1B
stackedBarData <- integratedData@meta.data %>% 
    group_by(sample,Cluster) %>% 
    summarize(count = n())


stackedBarPlot <- ggplot(stackedBarData, aes(fill=sample, y=count, x=as.factor(Cluster))) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_discrete(labels = c('Recurrent Glioblastoma 1 ','Oligodendroglioma 2 ','Glioblastoma 3','Recurrent Glioblastoma 4','Astrocytoma 5',
                                   'Astrocytoma 6','Diffuse Midline Glioma 7','Oligodendroglioma 8',
                                   'Glioblastoma 9','Glioblastoma 10','Oligodendroglioma 11')) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color='black'),
          strip.text = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.text.y = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          legend.position = 'none',
          # legend.text = element_text(size = 26, family = 'Arial Narrow', face = 'bold'),
          # legend.title = element_text(size = 30, family = 'Arial Narrow', face = 'bold'),
          plot.title = element_text(size = 30, family = 'Arial Narrow', face = 'bold', hjust = 0.5)) +
    xlab('Cluster') +
    labs(fill='Sample') +
    ggtitle('Cluster Distribution')
ggsave(stackedBarPlot, filename = 'integratedClusters_stackedBarPlot.png')


# Figure 1C
integratedData <- Seurat::FindNeighbors(integratedData)
integratedData <- Seurat::FindClusters(integratedData, resolution = 0.2)
integratedData@meta.data$Cluster <- as.numeric(integratedData@meta.data$integrated_snn_res.0.2)
clusterUMAP <- DimPlot(integratedData, reduction = 'umap', group.by = 'Cluster') +
    scale_color_manual(values=distinctColors) +
    ggtitle('Integrated Clusters') +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color='black'),
          strip.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.text.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.x = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'))
ggsave(clusterUMAP, filename = 'integratedClusters_UMAP_highlightedClusters.png')

# Figure 1D
clusterAssignments <- integratedData@meta.data %>% select(Cluster,sample)
clusterAssignments$Barcode <- row.names(clusterAssignments)
clusterAssignments$Barcode <- gsub('_.*','',clusterAssignments$Barcode)

sampleList <- unique(clusterAssignments$sample)
listOfClusters <- clusterAssignments %>% group_by(sample) %>% group_split()
listOfC <- list()
for(i in 1:length(listOfClusters)){
    df <- data.frame(listOfClusters[[i]])
    row.names(df) <- df$Barcode
    listOfC[[i]] <- df
}

listD <- list()
for(i in 1:length(listOfSamples2)){
    s <- listOfSamples2[[i]]
    c <- listOfC[[i]]
    s@meta.data <- transform(merge(s@meta.data,c,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
    listD[[i]] <- s
}

samples <- c('FFD1','FFSlideDB1', 'FFSlideDC1', 'FFSlideDD1', 'FFSlideEB1','FFSlideEC1','FFSlideED1','FFSlideEE1','FFSlideGB1','FFSlideGD1','FFSlideGE1')

labels <- c('Recurrent Glioblastoma 1 ','Oligodendroglioma 2 ','Glioblastoma 3','Recurrent Glioblastoma 4','Astrocytoma 5',
            'Astrocytoma 6','Diffuse Midline Glioma 7','Oligodendroglioma 8',
            'Glioblastoma 9','Glioblastoma 10','Oligodendroglioma 11')
listOfS <- list()
for(i in 1:length(listD)){
    s <- SpatialDimPlot(listD[[i]], group.by = 'Cluster',image.alpha = 0,stroke=0) +
        scale_fill_manual(values = distinctColors[1:9]) +
        # ggtitle(labels[[i]]) +
        theme(legend.key = element_rect(colour = NA, fill = NA),
              legend.box.background = element_blank(),
              legend.position='none',
              legend.text = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
              legend.title = element_text(size = 15, family = 'Arial Narrow', face = 'bold'),
              plot.title = element_text(size = 10, family = 'Arial Narrow', face = 'bold'),
              axis.line = element_line(color = 'black'), aspect.ratio = 1) +
        guides(fill = guide_legend(override.aes = list(size=5)))
    ggsave(s,filename = paste0('~/',samples[[i]],'_integratedClusters.png'))
}

# Figure 1E
Idents(integratedData) <- 'Cluster'
DefaultAssay(integratedData) <- 'Spatial'
integratedData <- NormalizeData(integratedData)
integratedData <- FindVariableFeatures(integratedData, selection.method = "vst", nfeatures = 2000)
integratedData <- ScaleData(integratedData, vars.to.regress = c('pctMito'))
markers <- FindAllMarkers(integratedData, assay = 'Spatial')
write.csv(markers,'integratedData_markerGenes.csv')

top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)


top10$cluster <- as.numeric(top10$cluster)

a <- DoHeatmap(integratedData,features = top10$gene, group.by = 'Cluster',
               angle = 0, group.colors = distinctColors[1:9], slot = 'scale.data') + NoLegend() +
    theme(legend.position='right',
          axis.text.y = element_text(size = 6, family = 'Arial Narrow', face = 'bold'),
          axis.title.y = element_text(size = 6, family = 'Arial Narrow', face = 'bold'),
          legend.text = element_text(size = 6, family = 'Arial Narrow', face = 'bold'),
          legend.title = element_text(size = 6, family = 'Arial Narrow', face = 'bold')) +
    guides(color = 'none')




# Figure 1F
integratedData@meta.data$type <- recode(integratedData@meta.data$sample, 
                                      "FFD1" = 'Recurrent Glioblastoma 1', 
                                      "FFSlideDB1" = 'Oligodendroglioma 2', 
                                      "FFSlideDC1" = 'Glioblastoma 3', 
                                      "FFSlideDD1" = 'Recurrent Glioblastoma 4',
                                      "FFSlideEB1" = 'Astrocytoma 5',
                                      "FFSlideEC1" = 'Astrocytoma 6',
                                      "FFSlideED1" = 'Diffuse Midline Glioma 7',
                                      "FFSlideEE1" = 'Oligodendroglioma 8',
                                      "FFSlideGB1" = 'Glioblastoma 9',
                                      "FFSlideGD1" = 'Glioblastoma 10',
                                      "FFSlideGE1" = 'Oligodendroglioma 11')

raviTranscriptionalModules <- read.csv('~/FourModules.tsv', sep='\t')
module1 <- raviTranscriptionalModules %>% filter(module == 'module_1') %>% pull(gene)
module2 <- raviTranscriptionalModules %>% filter(module == 'module_2') %>% pull(gene)
module3 <- raviTranscriptionalModules %>% filter(module == 'module_3') %>% pull(gene)
module4 <- raviTranscriptionalModules %>% filter(module == 'module_4') %>% pull(gene)
integratedData <- AddModuleScore(object = integratedData,
                                    features = list(module1,
                                                    module2,
                                                    module3,
                                                    module4),
                                    search = FALSE,
                                    name = c('Tumor Core', 'Vascular Niche', 'Invasive Niche', 'Hypoxic Niche'))

expressionTable <- as.matrix(integratedData@meta.data[c('Tumor.Core1','Vascular.Niche2','Invasive.Niche3','Hypoxic.Niche4')])
colnames(expressionTable) <- gsub('[0-9]','',colnames(expressionTable))

samples <- c("Recurrent Glioblastoma 1" = "#F8766D", "Oligodendroglioma 2" = "#DE8C00", "Glioblastoma 3" = "#B79F00",
             "Recurrent Glioblastoma 4" = "#7CAE00","Astrocytoma 5" = "#00BA38", "Astrocytoma 6" = "#00C08B",
             "Diffuse Midline Glioma 7" = "#00B4F0", "Oligodendroglioma 8" = "#619CFF", "Glioblastoma 9" = "#C77CFF",
             "Glioblastoma 10" = "#F564E3", "Oligodendroglioma 11" = "#FF64B0")

activity=list(Sample=samples)
row_ha <- rowAnnotation(Sample = integratedData@meta.data$type)
split = rep(1:4, each = 1)
ha = HeatmapAnnotation(
    empty = anno_empty(border = FALSE),
    foo = anno_block(gp = gpar(fill = c('gray93'), fontsize = 5), 
                     labels = c('Tumor Core', 'Vascular Niche', 'Invasive Niche', 'Hypoxic Niche'))
)

test <- as.matrix(integratedData@meta.data[c(0,11)])
colnames(test)[1] <- 'Sample'
row_ha = rowAnnotation(Sample = test, 
                       col = activity,
                       annotation_name_gp= gpar(fontsize = 10),
                       annotation_legend_param = list(Sample = list( at = unique(integratedData@meta.data$type))))
ht = Heatmap(expressionTable,
             name = 'Module Score',
             column_split = split,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             cluster_rows = FALSE,
             column_title = NULL,
             cluster_columns = FALSE,
             right_annotation = row_ha)

png(file="GBM_heatmap1.png", width = 650, height = 650)
pushViewport(viewport(gp = gpar(fontfamily = "Arial Narrow")))
draw(ht, padding = unit(c(2, 2, 2, 40), "mm"))
dev.off()


# Figure 1G
renRegionalPrograms <- read.csv('~/markerSet.tsv',sep='\t')
set1 <- renRegionalPrograms %>% filter(Regional.Programm == 'Radial.Glia') %>% pull(Gene)
set2 <- renRegionalPrograms %>% filter(Regional.Programm == 'Reactive.Immune') %>% pull(Gene)
set3 <- renRegionalPrograms %>% filter(Regional.Programm == 'Regional.NPC') %>% pull(Gene)
set4 <- renRegionalPrograms %>% filter(Regional.Programm == 'Regional.OPC') %>% pull(Gene)
set5 <- renRegionalPrograms %>% filter(Regional.Programm == 'Reactive.Hypoxia') %>% pull(Gene)

integratedData <- AddModuleScore(object = integratedData,
                                    features = list(set1,
                                                    set2,
                                                    set3,
                                                    set4,
                                                    set5),
                                    search = FALSE,
                                    name = c('Radial.Glia', 'Reactive.Immune', 'Regional.NPC', 'Regional.OPC','Reactive.Hypoxia'))

expressionTable <- as.matrix(integratedData@meta.data[c('Radial.Glia1', 'Reactive.Immune2', 'Regional.NPC3', 'Regional.OPC4','Reactive.Hypoxia5')])
colnames(expressionTable) <- gsub('[0-9]','',colnames(expressionTable))
activity=list(Sample=samples)
row_ha <- rowAnnotation(Sample = integratedData@meta.data$type)
split = rep(1:5, each = 1)
ha = HeatmapAnnotation(
    empty = anno_empty(border = FALSE),
    foo = anno_block(gp = gpar(fill = c('gray93'), fontsize = 5), 
                     labels = c('Radial Glia', 'Reactive Immune', 'Regional NPC', 'Regional OPC','Reactive Hypoxia'))
)

test <- as.matrix(integratedData@meta.data[c(0,11)])
colnames(test)[1] <- 'Sample'
row_ha = rowAnnotation(Sample = test, col = activity, 
                       annotation_legend_param = list(Sample = list( at = unique(integratedData@meta.data$type))))


ht = Heatmap(expressionTable,
             name = 'Module Score',
             column_split = split,
             top_annotation = ha,
             show_row_names = FALSE,
             show_column_names = FALSE,
             cluster_rows = FALSE,
             column_title = NULL,
             cluster_columns = FALSE,
             right_annotation = row_ha)

pushViewport(viewport(gp = gpar(fontfamily = "Arial Narrow")))
draw(ht, padding = unit(c(2, 2, 2, 40), "mm"))

png(file="GBM_regionalPrograms_heatmap2.png", width = 800, height = 800)
pushViewport(viewport(gp = gpar(fontfamily = "Arial Narrow")))
draw(ht, padding = unit(c(2, 2, 2, 40), "mm"))
dev.off()


