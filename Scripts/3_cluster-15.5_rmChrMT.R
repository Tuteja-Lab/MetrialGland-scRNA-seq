library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(1234)

load("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/2_mergedSamples/gd15.5.RNAmerged.rda")

data <- gd15.5.merged

#remove genes on chromosome MT
mart <- biomaRt::useMart("ensembl", host="http://sep2019.archive.ensembl.org", dataset = "rnorvegicus_gene_ensembl") #to make sure that we use the correct version of Ensembl
rnor6 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"), mart = mart) #build attribute table
mt.genes <- subset(rnor6, rnor6$chromosome_name == "MT")$external_gene_name


counts <- GetAssayData(data, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% mt.genes)),]
data <- subset(data, features = rownames(counts))


# Normalize data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(data)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/GD15.5/gd15.5_varFeatures.pdf", width = 10, height = 8)
#plot1 + plot2
#dev.off()

# Scale data
data <- ScaleData(data)

# PCA
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 100)

# Determining dimensions
#data <- JackStraw(data, num.replicate = 100, dims = 100)
#data <- ScoreJackStraw(data, dims  = 1:100)

#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_JackStrawPlot-dims100.pdf", width = 10, height = 8)
#JackStrawPlot(data, dims = 1:100)
#dev.off()

#pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_ElbowPlot.pdf", width = 10, height = 8)
#ElbowPlot(data, ndims = 100)
#dev.off()


# Cluster cells + UMAP
dim <- 72
data <- FindNeighbors(data, dims = 1:dim)


data <- FindClusters(data, resolution = 0.8)
data <- RunUMAP(data, dims = 1:dim)
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_dimPlot-res0.8.pdf", width = 10, height = 8)
DimPlot(data, reduction = "umap")
dev.off()
table(Idents(data))

save(data, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/3_cluster_DEG/rmChrMTgenes/GD15.5/gd15.5_res0.8.rda")







