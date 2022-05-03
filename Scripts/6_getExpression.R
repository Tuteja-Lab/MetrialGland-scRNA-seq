library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(1234)

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/4_twoTimepoints/rmChrMTgenes/DEG/twotimepoint.rda")
twotimepoint <- RenameIdents(twotimepoint, `23` = "gd19.5_23", `24` = "gd19.5_24", `25` = "gd19.5_25", `26` = "gd19.5_26")
data.expr <- AverageExpression(twotimepoint, slot  = "data")
expr <- data.expr$RNA[,order(colnames(data.expr$RNA))]
expr[1:5,1:5]

write.table(expr, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/4_twoTimepoints/rmChrMTgenes/twotimepointExpr_slotData.txt", sep = "\t", row.names = T, quote = F)

#get percent of cells that express a gene
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_res0.8.rda")
gd15.5 <- data
DefaultAssay(gd15.5) <- "RNA"
gd15.5 <- RenameIdents(gd15.5, `0` = "gd15.5_0", `1` = "gd15.5_1", `2` = "gd15.5_2", `3` = "gd15.5_3", `4` = "gd15.5_4", `5` = "gd15.5_5", `6` = "gd15.5_6-NKs", `7` = "gd15.5_7", `8` = "gd15.5_8", `9` = "gd15.5_9", `10` = "gd15.5_10", `11` = "gd15.5_11", `12` = "gd15.5_12", `13` = "gd15.5_13", `14` = "gd15.5_14", `15` = "gd15.5_15", `16` = "gd15.5_16", `17` = "gd15.5_17", `18` = "gd15.5_18", `19` = "gd15.5_19", `20` = "gd15.5_20", `21` = "gd15.5_21", `22` = "gd15.5_22", `23` = "gd15.5_23-TBC", `24` = "gd15.5_24", `25` = "gd15.5_25", `26` = "gd15.5_26", `27` = "gd15.5_27", `28` = "gd15.5_28")
gd15.5.pct <- DotPlot(gd15.5, features = rownames(gd15.5))
data15.5 <- subset(gd15.5.pct$data, gd15.5.pct$data$pct.exp >= 10)


load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_res0.8.rda")
gd19.5 <- data2
gd19.5 <- RenameIdents(gd19.5, `0` = "gd19.5_0", `1` = "gd19.5_1", `2` = "gd19.5_2", `3` = "gd19.5_3", `4` = "gd19.5_4", `5` = "gd19.5_5", `6` = "gd19.5_6-TBC", `7` = "gd19.5_7", `8` = "gd19.5_8", `9` = "gd19.5_9-NKs", `10` = "gd19.5_10", `11` = "gd19.5_11", `12` = "gd19.5_12", `13` = "gd19.5_13", `14` = "gd19.5_14", `15` = "gd19.5_15", `16` = "gd19.5_16", `17` = "gd19.5_17", `18` = "gd19.5_18", `19` = "gd19.5_19", `20` = "gd19.5_20", `21` = "gd19.5_21", `22` = "gd19.5_22", `23` = "gd19.5_23", `24` = "gd19.5_24", `25` = "gd19.5_25", `26` = "gd19.5_26")

gd19.5.pct <- DotPlot(gd19.5, features = rownames(gd19.5@assays$RNA), assay = "RNA")
data19.5 <- subset(gd19.5.pct$data, gd19.5.pct$data$pct.exp >= 10)

expr <- read.table("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/4_twoTimepoints/rmChrMTgenes/twotimepointExpr_slotData.txt", sep = "\t", header = T)

expr1 <- subset(expr, rownames(expr) %in% data15.5$features.plot & rownames(expr) %in% data19.5$features.plot)

write.table(expr1, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/4_twoTimepoints/rmChrMTgenes/twotimepointExpr_10pctBoth_slotData.txt", sep = "\t", row.names = T, quote = F)

expr1 <- subset(expr, rownames(expr) %in% data19.5$features.plot)

write.table(expr1, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/4_twoTimepoints/rmChrMTgenes/twotimepointExpr_10pctGD19.5_slotData.txt", sep = "\t", row.names = T, quote = F)






