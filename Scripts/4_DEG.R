library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(1234)



#load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD15.5/gd15.5_res0.5.rda")

#gd15.5.markers.res0.5.wilcoxin <- FindAllMarkers(data, only.pos = TRUE, logfc.threshold=0)
#gd15.5.markers.res0.5.mast <- FindAllMarkers(data, only.pos = TRUE, logfc.threshold=0, test.use = "MAST")

#save(gd15.5.markers.res0.5.wilcoxin, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD15.5/gd15.5.markers.res0.5.wilcoxin.rda")
#save(gd15.5.markers.res0.5.mast, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD15.5/gd15.5.markers.res0.5.mast.rda")




load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_res0.8.rda")

gd15.5.markers.res0.8.wilcoxin <- FindAllMarkers(data, only.pos = TRUE, logfc.threshold=0)
#gd15.5.markers.res0.8.mast <- FindAllMarkers(data, only.pos = TRUE, logfc.threshold=0, test.use = "MAST")

save(gd15.5.markers.res0.8.wilcoxin, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5.markers.res0.8.wilcoxin.rda")
#save(gd15.5.markers.res0.8.mast, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD15.5/gd15.5.markers.res0.8.mast.rda")







