library(Seurat)
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD15.5_1-filtered.rda")
gd15.5_1 <- rats.rna
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD15.5_2-filtered.rda")
gd15.5_2 <- rats.rna
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD15.5_3-filtered.rda")
gd15.5_3 <- rats.rna
gd15.5.merged <- merge(gd15.5_1, y = c(gd15.5_2, gd15.5_3), add.cell.ids = c("15.5_1", "15.5_2", "15.5_3"), project = "RNA_15.5")
gd15.5.merged
save(gd15.5.merged, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/2_mergedSamples/gd15.5.RNAmerged.rda")




load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_4-filtered.rda")
gd19.5_4 <- rats.rna
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_5-filtered.rda")
gd19.5_5 <- rats.rna
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_6-filtered.rda")
gd19.5_6 <- rats.rna
gd19.5.merged <- merge(gd19.5_4, y = c(gd19.5_5, gd19.5_6), add.cell.ids = c("19.5_4", "19.5_5", "19.5_6"), project = "RNA_19.5")
save(gd19.5.merged, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/2_mergedSamples/gd19.5.RNAmerged.rda")
