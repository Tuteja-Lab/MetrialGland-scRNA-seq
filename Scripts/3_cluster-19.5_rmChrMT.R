library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(1234)

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/2_mergedSamples/integrationGD19.5/rmChrMTgenes/gd19.5-4-5-6-7.integrated.rda")

data <- gd19.5.combined

DefaultAssay(data) <- "integrated"

#remove genes on chromosome MT
#mart <- biomaRt::useMart("ensembl", host="http://sep2019.archive.ensembl.org", dataset = "rnorvegicus_gene_ensembl") #to make sure that we use the correct version of Ensembl
#rnor6 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"), mart = mart) #build attribute table
#mt.genes <- subset(rnor6, rnor6$chromosome_name == "MT")$external_gene_name
#existingMTgenes <- intersect(rownames(data), mt.genes)

#counts <- GetAssayData(data, assay = "RNA")
#counts <- counts[-(which(rownames(counts) %in% existingMTgenes)),]
#data <- subset(data, features = rownames(counts))



# Normalize data
#data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
#data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(data)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD19.5/gd19.5_varFeatures.pdf", width = 10, height = 8)
#plot1 + plot2
#dev.off()

# Scale data
data <- ScaleData(data)

# PCA
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 100)

# Determining dimensions
#data <- JackStraw(data, num.replicate = 100, dims = 100)
#data <- ScoreJackStraw(data, dims  = 1:100)

#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_JackStrawPlot-dims100.pdf", width = 10, height = 8)
#JackStrawPlot(data, dims = 1:100)
#dev.off()

#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_ElbowPlot.pdf", width = 10, height = 8)
#ElbowPlot(data, ndims = 100)
#dev.off()


# Cluster cells + UMAP
dim <- 77
data <- FindNeighbors(data, dims = 1:dim)

#test resolutions
data2 <- FindClusters(data, resolution = 0.8)
data2 <- RunUMAP(data2, dims = 1:dim)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_dimPlot-res0.8.pdf", width = 10, height = 8)
DimPlot(data2, reduction = "umap")
dev.off()
table(Idents(data2))

save(data2, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD19.5/gd19.5-4-5-6-7_res0.8.rda")







