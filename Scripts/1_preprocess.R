library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#load dataset
rna.data <- Read10X(data.dir = "/work/LAS/geetu-lab/hhvu/project3_scATAC/data/scRNA-seq/GD19.5-RNA/RNA-19.5-6")

# Initialize the Seurat object with the raw (non-normalized data)
rats.rna <- CreateSeuratObject(counts = rna.data, project = "RNA", min.cells = 1)
rats.rna

# The [[ operator can add columns to object metadata. This is a great place to stash QC stat
mart <- biomaRt::useMart("ensembl", host="http://sep2019.archive.ensembl.org", dataset = "rnorvegicus_gene_ensembl") #to make sure that we use the correct version of Ensembl
rnor6 <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"), mart = mart) #build attribute table
mt.genes <- subset(rnor6, rnor6$chromosome_name == "MT")$external_gene_name
existingMTgenes <- intersect(rownames(rats.rna), mt.genes)


rats.rna[["percent.mt"]] <- PercentageFeatureSet(rats.rna, features = existingMTgenes)


# Visualize QC metrics as a violin plot
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_QC/GD19.5_6_comboPlots.pdf", width = 10, height = 8)
VlnPlot(rats.rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

print("nFeature_RNA percentile")
quantile(rats.rna[["nFeature_RNA"]][, 1], probs = seq(0, 1, 0.25))

print("nCount_RNA percentile")
quantile(rats.rna[["nCount_RNA"]][, 1], probs = seq(0, 1, 0.25))

print("percent.mt percentile")
quantile(rats.rna[["percent.mt"]][, 1], probs = seq(0, 1, 0.25))

plot2 <- FeatureScatter(rats.rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(rats.rna, feature1 = "nFeature_RNA", feature2 = "percent.mt")
pdf("/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_QC/GD19.5_6_featureScatter.pdf", width = 10, height = 8)
plot2 + plot3
dev.off()



#filtering
#test1 <- subset(rats.rna, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 15)
#test1

#test1 <- subset(rats.rna, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 20)
#test1



rats.rna <- subset(rats.rna, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 20)
#rats.rna

save(rats.rna, file = "/work/LAS/geetu-lab/hhvu/project3_scATAC/scRNA-seq-analysis/1_preprocess/GD19.5_6-filtered.rda")
