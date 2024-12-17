library(Seurat) # version 5.1.0
library(dplyr) # version 1.1.4  

setwd('~Figures')

GlioblastomaA1 <- Load10X_Spatial('~/GlioblastomaA1/outs/')
GlioblastomaA1 <- NormalizeData(GlioblastomaA1)

a1 <- SpatialFeaturePlot(GlioblastomaA1, features = 'VEGFA',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a1, file = 'SupplementaryFigure6_a_VEGFA.png', height = 7, width = 7)
a2 <- SpatialFeaturePlot(GlioblastomaA1, features = 'HIF1A',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a2, file = 'SupplementaryFigure6_b_HIF1A.png', height = 7, width = 7)
a3 <- SpatialFeaturePlot(GlioblastomaA1, features = 'PDGFRA',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a3, file = 'SupplementaryFigure6_c_PDGFRA.png', height = 7, width = 7)
a4 <- SpatialFeaturePlot(GlioblastomaA1, features = 'NDRG1',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a4, file = 'SupplementaryFigure6_d_NDRG1.png', height = 7, width = 7)
a5 <- SpatialFeaturePlot(GlioblastomaA1, features = 'DDIT4',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a5, file = 'SupplementaryFigure6_e_DDIT4.png', height = 7, width = 7)
a6 <- SpatialFeaturePlot(GlioblastomaA1, features = 'CD44',slot = 'data',image.alpha=0,pt.size.factor = 2.5) +
    theme(legend.position='right',
          legend.title = element_text(face='bold.italic',family='Arial Narrow',size=15),
          legend.text = element_text(face='bold',family='Arial Narrow',size=15))
ggsave(a6, file = 'SupplementaryFigure6_f_CD44.png', height = 7, width = 7)

