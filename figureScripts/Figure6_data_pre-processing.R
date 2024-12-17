library(Seurat)

library(SeuratData)

library(Matrix)

library(patchwork)

library(dplyr)

library(loupeR)

library(rhdf5)

library(ggplot2)

library(harmony)

options(future.globals.maxSize = 2048 * 1024^2)  # Set limit to 2 GB

sample="spa5"

expression_data <- Read10X(data.dir = "Craig_SPA5_D/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA5_D/outs/"

visium_data <- CreateSeuratObject(counts = expression_data)

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

sample="spa1"

expression_data <- Read10X(data.dir = "Craig_SPA1_D/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA1_D/outs/"

visium_data <- CreateSeuratObject(counts = expression_data)

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

sample="spa3"

expression_data <- Read10X(data.dir = "Craig_SPA3_D/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA3_D/outs/"

visium_data <- CreateSeuratObject(counts = expression_data)

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

sample="spa8"

expression_data <- Read10X(data.dir = "Craig_SPA8_A/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA8_A/outs/"

visium_data <- CreateSeuratObject(counts = expression_data)

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

sample="spa2"

expression_data <- Read10X(data.dir = "Craig_SPA2_A/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA2_A/outs/"

visium_data <- CreateSeuratObject(counts = expression_data)

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

sample="spa7"

expression_data <- Read10X(data.dir = "Craig_SPA7_D/outs/filtered_feature_bc_matrix/")

img_data_dir <- "Craig_SPA7_D/outs/"

spatial_img <- Read10X_Image(image.dir = file.path(img_data_dir, "spatial"))

visium_data <- CreateSeuratObject(counts = expression_data, assay = "Spatial")

visium_data[["slice1"]] <- spatial_img



visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE,)

visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)

visium_data <- FindClusters(visium_data, resolution=0.5, verbose = FALSE)

visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)

p<-p1+p2

ggsave(paste0(sample,".umap.png"), plot = p, width = 6, height = 4, units = "in", dpi = 300)

cluster_markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- visium_data@assays$SCT@counts

clusters <- select_clusters(visium_data)

projections <- select_projections(visium_data)

write.csv(cluster_markers, paste0(sample,".seurat.clusterMarkers.seurat.tsv"))

write.csv(clusters, paste0(sample,".clusters.seurat.tsv"))

saveRDS(visium_data,file = paste0(sample,".seurat.RDS"))

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0(sample,".seurat.cloupe"),force = TRUE)

options(future.globals.maxSize = 2048 * 1024^2)  # Set limit to 2 GB



spa1 <- readRDS("spa1.seurat.RDS")

spa2 <- readRDS("spa2.seurat.RDS")

spa3 <- readRDS("spa3.seurat.RDS")

spa5 <- readRDS("spa5.seurat.RDS")

spa7 <- readRDS("spa7.seurat.RDS")

spa8 <- readRDS("spa8.seurat.RDS")

spa1$orig.ident <- "spa1"

spa2$orig.ident <- "spa2"

spa3$orig.ident <- "spa3"

spa5$orig.ident <- "spa5"

spa7$orig.ident <- "spa7"

spa8$orig.ident <- "spa8"

spa1 <- SCTransform(spa1, assay = "Spatial",verbose = FALSE,variable.features.n = 10000)

spa2 <- SCTransform(spa2, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

spa3 <- SCTransform(spa3, assay = "Spatial",verbose = FALSE,variable.features.n = 10000)

spa5 <- SCTransform(spa5, assay = "Spatial", verbose = FALSE,variable.features.n = 10000)

spa7 <- SCTransform(spa7, assay = "Spatial",verbose = FALSE,variable.features.n = 10000)

spa8 <- SCTransform(spa8,assay = "Spatial",verbose = FALSE,variable.features.n = 10000)

spa1 <- FindVariableFeatures(spa1, selection.method = "vst", nfeatures = 4000)

spa2 <- FindVariableFeatures(spa2, selection.method = "vst", nfeatures = 4000)

spa3 <- FindVariableFeatures(spa3, selection.method = "vst", nfeatures = 4000)

spa5 <- FindVariableFeatures(spa5, selection.method = "vst", nfeatures = 4000)

spa7 <- FindVariableFeatures(spa7, selection.method = "vst", nfeatures = 4000)

spa8 <- FindVariableFeatures(spa8, selection.method = "vst", nfeatures = 4000)

var_features_spa1 <- VariableFeatures(spa1)

var_features_spa2 <- VariableFeatures(spa2)

var_features_spa3 <- VariableFeatures(spa3)

var_features_spa5 <- VariableFeatures(spa5)

var_features_spa7 <- VariableFeatures(spa7)

var_features_spa8 <- VariableFeatures(spa8)

combined <- merge(spa1, y = c(spa2, spa3, spa5,spa7,spa8), add.cell.ids = c("spa1", "spa2", "spa3","spa5","spa7","spa8"))

DoHeatmap(combined, features = c("ADM","BCAN","NCAN"), size = 3) + NoLegend()

combined$orig.ident <- as.factor(combined$orig.ident)

combined_var_features <- unique(c(var_features_spa1, var_features_spa2, var_features_spa3, var_features_spa5,var_features_spa7,var_features_spa8))

VariableFeatures(combined) <- combined_var_features

combined <- RunPCA(combined,  npcs = 30, verbose = FALSE)

combined$orig.ident <- as.factor(combined$orig.ident)

combined <- RunHarmony(combined, group.by.vars = "orig.ident")

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)

combined <- FindClusters(combined, resolution = .19)

DimPlot(combined, reduction = "umap", group.by = "orig.ident")

DimPlot(combined, reduction = "umap", label = TRUE)

umap_coordinates <- Embeddings(combined, "umap")



# Extract cluster labels

cluster_labels <- Idents(combined)

umap_df <- data.frame(umap_coordinates, cluster = cluster_labels)

write.csv(umap_df, file = "umap_coordinates_clusters.csv", row.names = TRUE)

# Extract the raw counts matrix

counts_matrix <- GetAssayData(combined, slot = "counts")

annotation <- data.frame(cell = colnames(counts_matrix),
                         
                         cell_type = Idents(combined))

p1 <- DimPlot(combined, reduction = "umap", label = TRUE) + NoLegend()

p2 <- SpatialDimPlot(combined, label = TRUE, label.size = 3)

p<-p1+p2

ggsave("combined.umap.png", plot = p, width = 6, height = 4, units = "in", dpi = 300)

combined <- PrepSCTFindMarkers(combined)



cluster_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

count_mat <- combined@assays$SCT@counts

clusters <- select_clusters(combined)

projections <- select_projections(combined)

write.csv(cluster_markers, "combined.clusterMarkers.csv")

write.csv(clusters, "combined.clusters.csv")

write.table(annotation, file="combined.annotations.tsv", sep="\t", quote=FALSE, row.names=FALSE)



saveRDS(combined,file = "combined.seurat.RDS")

create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = "combined.harmony.cloupe",force = TRUE)


saveRDS(counts_matrix,file = "combined.counts_matrix.RDS")

saveRDS(combined,file = "combined.RDS")

combined <- readRDS("combined.seurat.RDS")

cluster_markers<-read.table(file="fig.combined.clusterMarkers.csv", sep = ",", header = TRUE)

top_markers <- cluster_markers %>%
    
    group_by(cluster) %>%
    
    top_n(n = 10, wt = avg_log2FC)



DoHeatmap(combined, features = top_markers$gene, size = 3) + NoLegend()

barcode_cluster_df <- read.csv("test.csv")

rownames(barcode_cluster_df) <- barcode_cluster_df$barcode

barcode_cluster_df$barcode <- NULL

combined <- AddMetaData(combined, metadata = barcode_cluster_df)
