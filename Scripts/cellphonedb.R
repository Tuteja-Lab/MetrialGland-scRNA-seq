library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
set.seed(1234)

# module load miniconda3/4.3.30-qdauveb
# conda create -n cpdb python=3.7
# module load gcc/7.3.0-xegsmw4
# module load r/4.0.2-py3-icvulwq
# pip install cellphonedb

# module load miniconda3/4.3.30-qdauveb
# module load gcc/7.3.0-xegsmw4
# module load r/4.0.2-py3-icvulwq
# source activate cpdb
# source deactivate cpdb

#gd15.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_res0.8.rda")
gd15.5 <- data
gd15.5$sample <- substr(names(gd15.5$orig.ident), 1, 6)
DefaultAssay(gd15.5) <- "RNA"
gd15.5 <- RenameIdents(gd15.5, `0` = "Other cells", `1` = "Other cells",
                       `2` = "Other cells", `3` = "Other cells", `4` = "Other cells",
                       `5` = "Other cells", `6` = "Natural killer cells",
                       `7` = "Macrophage cells", `8` = "Macrophage cells",
                       `9` = "Natural killer cells", `10` = "Endothelial cells",
                       `11` = "Macrophage cells", `12` = "Other cells",
                       `13` = "Endothelial cells", `14` = "Smooth muscle cells",
                       `15` = "Endothelial cells", `16` = "Other cells",
                       `17` = "Smooth muscle cells", `18` = "Other cells",
                       `19` = "Natural killer cells", `20` = "Other cells",
                       `21` = "Other cells", `22` = "Other cells",
                       `23` = "Invasive trophoblast cells", `24` = "Other cells",
                       `25` = "Other cells", `26` = "Smooth muscle cells",
                       `27` = "Other cells", `28` = "Other cells")
Idents(gd15.5) <- factor(Idents(gd15.5), levels = c("Other cells", "Macrophage cells", "Smooth muscle cells",
                                                   "Endothelial cells", "Natural killer cells", "Invasive trophoblast cells"))
load(file = "/work/LAS/geetu-lab/hhvu/combine-test-expression1.Rdata")
ratHumanOrthologs <- dataset$GRCH38$ratHumanOrthologs
humanGeneMapping <- dataset$GRCH38$humanGeneMapping

genes <- rownames(gd15.5@assays$RNA@data)
genes <- toupper(genes)
ratHumanOrthologs <- ratHumanOrthologs[ratHumanOrthologs$Gene.name %in% genes,]

d <- gd15.5@assays$RNA@data
d <- d[toupper(rownames(d)) %in% ratHumanOrthologs$Gene.name,]
# Save normalised counts
#writeMM(d, file = '/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd15.5matrix.mtx')

# save gene and cell names
d <- as.data.frame(d)
d$genes <- toupper(rownames(d))
d <- inner_join(d, ratHumanOrthologs[,c("Gene.name", "Human.gene.name")], by = c("genes" = "Gene.name"))
rownames(d) <- d$Human.gene.name
d <- d[,1:(ncol(d)-2)]
write(x = rownames(d), file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd15.5features.tsv")
write(x = colnames(d), file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd15.5barcodes.tsv")
write.table(d, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd15.5counts.txt", sep = '\t', quote = F)

gd15.5@meta.data$idents <- Idents(gd15.5)
table(gd15.5@meta.data$idents)
gd15.5@meta.data$Cell = rownames(gd15.5@meta.data)
df <- gd15.5@meta.data[, c('Cell', 'idents')]
write.table(df, file ='/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd15.5_meta.tsv', sep = '\t', quote = F, row.names = F)


#gd19.5
load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_res0.8.rda")
gd19.5 <- data2
gd19.5 <- RenameIdents(gd19.5, `0` = "Other cells", `1` = "Macrophage cells",
                       `2` = "Other cells", `3` = "Other cells", `4` = "Other cells",
                       `5` = "Other cells", `6` = "Invasive trophoblast cells",
                       `7` = "Macrophage cells", `8` = "Smooth muscle cells",
                       `9` = "Natural killer cells", `10` = "Endothelial cells",
                       `11` = "Macrophage cells", `12` = "Other cells",
                       `13` = "Other cells", `14` = "Macrophage cells",
                       `15` = "Other cells", `16` = "Other cells",
                       `17` = "Macrophage cells", `18` = "Macrophage cells",
                       `19` = "Other cells", `20` = "Endothelial cells",
                       `21` = "Other cells", `22` = "Other cells",
                       `23` = "Other cells", `24` = "Other cells",
                       `25` = "Other cells", `26` = "Smooth muscle cells")
Idents(gd19.5) <- factor(Idents(gd19.5), levels = c("Other cells", "Macrophage cells", "Smooth muscle cells",
                                                   "Endothelial cells", "Natural killer cells", "Invasive trophoblast cells"))
DefaultAssay(gd19.5) <- "RNA"
load(file = "/work/LAS/geetu-lab/hhvu/combine-test-expression1.Rdata")
ratHumanOrthologs <- dataset$GRCH38$ratHumanOrthologs
humanGeneMapping <- dataset$GRCH38$humanGeneMapping

genes <- rownames(gd19.5@assays$RNA@data)
genes <- toupper(genes)
ratHumanOrthologs <- ratHumanOrthologs[ratHumanOrthologs$Gene.name %in% genes,]

d <- gd19.5@assays$RNA@data
d <- d[toupper(rownames(d)) %in% ratHumanOrthologs$Gene.name,]

# save gene and cell names
d <- as.data.frame(d)
d$genes <- toupper(rownames(d))
d <- inner_join(d, ratHumanOrthologs[,c("Gene.name", "Human.gene.name")], by = c("genes" = "Gene.name"))
rownames(d) <- d$Human.gene.name
d <- d[,1:(ncol(d)-2)]
write(x = rownames(d), file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd19.5features.tsv")
write(x = colnames(d), file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd19.5barcodes.tsv")
write.table(d, "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd19.5counts.txt", sep = '\t', quote = F)

gd19.5@meta.data$idents <- Idents(gd19.5)
table(gd19.5@meta.data$idents)
gd19.5@meta.data$Cell = rownames(gd19.5@meta.data)
df <- gd19.5@meta.data[, c('Cell', 'idents')]
write.table(df, file ='/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/8_cellphonedb/gd19.5_meta.tsv', sep = '\t', quote = F, row.names = F)
