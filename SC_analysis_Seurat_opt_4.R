library(tidyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library("tidyverse")
library(reshape2)
library(ggpubr)
library(DESeq2)
library(edgeR)
library(ggfortify)
library(factoextra)
library(tximport)
library(gplots)
library(tools)
library(data.table)
library(purrr)
library(eply)
library(tibble)


# Load the Tabula Muris spleen dataset

spl_m_data<- Read10X(data.dir = "Spleen-10X_P4_7/") #data folder

# Create Seurat object

spl_m <- CreateSeuratObject(counts = spl_m_data, project = "TM_Spleen", min.cells = 3, min.features = 200)
spl_m
spl_m<- spl_m[, sample(colnames(spl_m), size =1000, replace=F)]

rm(spl_m_data)

# QC and selecting cells for further analysis

spl_m[["percent.mt"]] <- PercentageFeatureSet(spl_m, pattern = "^MT-")

VlnPlot(spl_m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(spl_m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


spl_m <- subset(spl_m, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)

FeatureScatter(spl_m, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalizing the data

spl_m <- NormalizeData(spl_m)

# Identification of highly variable features

spl_m <- FindVariableFeatures(spl_m, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(spl_m), 10)

plot1 <- VariableFeaturePlot(spl_m)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

# Scaling the data

all.genes <- rownames(spl_m)
spl_m <- ScaleData(spl_m, features = all.genes)

# Run PCA

spl_m <- RunPCA(spl_m, features = VariableFeatures(object = spl_m))

print(spl_m[["pca"]], dims = 1:3, nfeatures = 3)
VizDimLoadings(spl_m, dims = 1:3, reduction = "pca")

DimPlot(spl_m, reduction = "pca")

# Determine the 'dimensionality' of the dataset

spl_m <- JackStraw(spl_m, num.replicate = 50)
spl_m <- ScoreJackStraw(spl_m, dims = 1:20)
#JackStrawPlot(spl_m, dims = 1:15)

ElbowPlot(spl_m)

# Cell clustering

spl_m <- FindNeighbors(spl_m, dims = 1:9)
spl_m <- FindClusters(spl_m, resolution = 0.4)

# Run Umap and tSNE

spl_m <- RunUMAP(spl_m, dims = 1:9)
DimPlot(spl_m, reduction = "umap", label='TRUE')


spl_m <- RunTSNE(spl_m, dims = 1:9)
DimPlot(spl_m, reduction = "tsne", label='TRUE')


# Finding differentially expressed features (cluster biomarkers)

spl_m.markers <- FindAllMarkers(spl_m, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers<-spl_m.markers %>% group_by(cluster) %>% top_n(n = 3)
gen <- top_markers$gene
gen
VlnPlot(spl_m, features = c("Ccl5","Lyz1","Lyz2","Gzma","Ifitm3","Gstm1","Cst3","S100a4","S100a6","Fcer1g", "Vpreb3"))
FeaturePlot(spl_m, features = c("Ccl5","Lyz1","Lyz2","Gzma","Ifitm3","Gstm1","Cst3","S100a4","S100a6","Fcer1g"))
top10

VlnPlot(spl_m, features = c("Zfp36l2","Zfp36l1","Zfp706","Unc93b1","Ublcp1","Ube2r2","Zap70","Znrf1","Zyx","Zfp36","Znhit1","Zfp593","Zmiz2", "Zfp330","Zmiz1","Zfp36l2"))  
VlnPlot(spl_m, features = c("Vpreb3","Zcchc11","Zfp36l1","Tsc22d3","Tnfrsf13c","Sypl","Zap70","Znrf1","Zdhhc20","Zyx","Zmiz2","Znhit1","Znhit1","Zmat2","Zranb2","Ywhaq"))

new.cluster.ids <- c("B cell","B_cell", "T cell", "monocytes", "NK cell","na")

names(new.cluster.ids) <- levels(spl_m)
spl_m <- RenameIdents(spl_m, new.cluster.ids)
DimPlot(spl_m, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

