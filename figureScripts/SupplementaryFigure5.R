library(dplyr) # version 1.1.4  
library(purrr) #  version 1.0.2
library(Seurat) # version 5.1.0


genes <- c('VEGFA','PDGFRA','HIF1A')
samples <- c('Glioblastoma A4','Glioblastoma A5','Oligodendroglioma A2')


GlioblastomaA3 <- Load10X_Spatial('~/Google Drive/Shared drives/10X_Alpha_Test_Phase_1_and_2/manuscript_revision3/spatialData/GlioblastomaA3 (1)/outs/')
GlioblastomaA5 <- Load10X_Spatial('~/Google Drive/Shared drives/10X_Alpha_Test_Phase_1_and_2/manuscript_revision3/spatialData/GlioblastomaA5/outs/')
OligodendrogliomaA3 <- Load10X_Spatial('~/Google Drive/Shared drives/10X_Alpha_Test_Phase_1_and_2/manuscript_revision3/spatialData/OligodendrogliomaA3/outs/')

GlioblastomaA3 <- NormalizeData(GlioblastomaA3)
SpatialFeaturePlot(GlioblastomaA3,features='VEGFA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(GlioblastomaA3,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(GlioblastomaA3,features='HIF1A',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')

GlioblastomaA5 <- NormalizeData(GlioblastomaA5)
SpatialFeaturePlot(GlioblastomaA5,features='VEGFA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(GlioblastomaA5,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(GlioblastomaA5,features='HIF1A',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')

OligodendrogliomaA3 <- NormalizeData(OligodendrogliomaA3)
SpatialFeaturePlot(OligodendrogliomaA3,features='VEGFA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(OligodendrogliomaA3,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')
SpatialFeaturePlot(OligodendrogliomaA3,features='HIF1A',image.alpha = 0, pt.size.factor = 2.3) + theme(legend.position = 'right')




























