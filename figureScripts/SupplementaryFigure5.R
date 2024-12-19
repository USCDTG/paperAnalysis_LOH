library(dplyr) # version 1.1.4  
library(purrr) #  version 1.0.2
library(Seurat) # version 5.1.0
library(viridis) # version 0.6.5
library(ggplot2) # version 3.5.1

GlioblastomaA3 <- Load10X_Spatial('~/GlioblastomaA3/outs/')
GlioblastomaA5 <- Load10X_Spatial('~/GlioblastomaA5/outs/')
OligodendrogliomaA3 <- Load10X_Spatial('~/OligodendrogliomaA3/outs/')

GlioblastomaA3 <- SCTransform(GlioblastomaA3, assay = "Spatial", verbose = FALSE)
a <- SpatialFeaturePlot(GlioblastomaA3,features='VEGFA',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(a, filename = 'SupplementaryFigure5_GlioblastomaA3_VEGFA.png',height=7,width=7)

b <- SpatialFeaturePlot(GlioblastomaA3,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(b, filename = 'SupplementaryFigure5_GlioblastomaA3_PDGFRA.png',height=7,width=7)

c <- SpatialFeaturePlot(GlioblastomaA3,features='HIF1A',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(c, filename = 'SupplementaryFigure5_GlioblastomaA3_HIF1A.png',height=7,width=7)

GlioblastomaA5 <- SCTransform(GlioblastomaA5, assay = "Spatial", verbose = FALSE)
d <- SpatialFeaturePlot(GlioblastomaA5,features='VEGFA',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(d, filename = 'SupplementaryFigure5_GlioblastomaA5_VEGFA.png',height=7,width=7)

e <- SpatialFeaturePlot(GlioblastomaA5,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(e, filename = 'SupplementaryFigure5_GlioblastomaA5_PDGFRA.png',height=7,width=7)

f <- SpatialFeaturePlot(GlioblastomaA5,features='HIF1A',image.alpha = 0, pt.size.factor = 2.6) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(f, filename = 'SupplementaryFigure5_GlioblastomaA5_HIF1A.png',height=7,width=7)

OligodendrogliomaA3 <- SCTransform(OligodendrogliomaA3, assay = "Spatial", verbose = FALSE)
g <- SpatialFeaturePlot(OligodendrogliomaA3,features='VEGFA',image.alpha = 0, pt.size.factor = 2.4) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(g, filename = 'SupplementaryFigure5_OligodendrogliomaA3_VEGFA.png',height=7,width=7)

h <- SpatialFeaturePlot(OligodendrogliomaA3,features='PDGFRA',image.alpha = 0, pt.size.factor = 2.4) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(h, filename = 'SupplementaryFigure5_OligodendrogliomaA3_PDGFRA.png',height=7,width=7)

i <- SpatialFeaturePlot(OligodendrogliomaA3,features='HIF1A',image.alpha = 0, pt.size.factor = 2.4) + 
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15)) +
    scale_fill_viridis(option='magma')

ggsave(i, filename = 'SupplementaryFigure5_OligodendrogliomaA3_HIF1A.png',height=7,width=7)




























