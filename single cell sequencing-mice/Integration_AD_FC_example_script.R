# example file

# Library usage
library(Seurat)
library(ggplot2)

setwd("C:/Users/Orchi/Box/DI/Single_cell/AD/GSE163577_RAW")

Sys.setenv('R_MAX_VSIZE'=32000000000)
memory.limit(240000000)

# 1. Data loading
HC1 <- readRDS("HC1.rds")
HC2 <- readRDS("HC2.rds")
HC3 <- readRDS("HC3.rds")
HC4 <- readRDS("HC4.rds")
HC5 <- readRDS("HC5.rds")
HC6 <- readRDS("HC6.rds")
HC7 <- readRDS("HC7.rds")
HC8 <- readRDS("HC8.rds")


HAD1 <- readRDS("HAD1.rds")
HAD2 <- readRDS("HAD2.rds")
HAD3 <- readRDS("HAD3.rds")
HAD4 <- readRDS("HAD4.rds")
HAD5 <- readRDS("HAD5.rds")
HAD6 <- readRDS("HAD6.rds")
HAD7 <- readRDS("HAD7.rds")
HAD8 <- readRDS("HAD8.rds")
HAD9 <- readRDS("HAD9.rds")

CC1 <- readRDS("CC1.rds")
CC2 <- readRDS("CC2.rds")
CC3 <- readRDS("CC3.rds")
CC4 <- readRDS("CC4.rds")

CAD1 <- readRDS("CAD1.rds")
CAD2 <- readRDS("CAD2.rds")
CAD3 <- readRDS("CAD3.rds")
CAD4 <- readRDS("CAD4.rds")

#mito
HC1[["percent.mt"]] <- PercentageFeatureSet(HC1, pattern = "^MT-")
HC2[["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^MT-")
HC3[["percent.mt"]] <- PercentageFeatureSet(HC3, pattern = "^MT-")
HC4[["percent.mt"]] <- PercentageFeatureSet(HC4, pattern = "^MT-")
HC5[["percent.mt"]] <- PercentageFeatureSet(HC5, pattern = "^MT-")
HC6[["percent.mt"]] <- PercentageFeatureSet(HC6, pattern = "^MT-")
HC7[["percent.mt"]] <- PercentageFeatureSet(HC7, pattern = "^MT-")
HC8[["percent.mt"]] <- PercentageFeatureSet(HC8, pattern = "^MT-")

HAD1[["percent.mt"]] <- PercentageFeatureSet(HAD1, pattern = "^MT-")
HAD2[["percent.mt"]] <- PercentageFeatureSet(HAD2, pattern = "^MT-")
HAD3[["percent.mt"]] <- PercentageFeatureSet(HAD3, pattern = "^MT-")
HAD4[["percent.mt"]] <- PercentageFeatureSet(HAD4, pattern = "^MT-")
HAD5[["percent.mt"]] <- PercentageFeatureSet(HAD5, pattern = "^MT-")
HAD6[["percent.mt"]] <- PercentageFeatureSet(HAD6, pattern = "^MT-")
HAD7[["percent.mt"]] <- PercentageFeatureSet(HAD7, pattern = "^MT-")
HAD8[["percent.mt"]] <- PercentageFeatureSet(HAD8, pattern = "^MT-")
HAD9[["percent.mt"]] <- PercentageFeatureSet(HAD9, pattern = "^MT-")

CC1[["percent.mt"]] <- PercentageFeatureSet(CC1, pattern = "^MT-")
CC2[["percent.mt"]] <- PercentageFeatureSet(CC2, pattern = "^MT-")
CC3[["percent.mt"]] <- PercentageFeatureSet(CC3, pattern = "^MT-")
CC4[["percent.mt"]] <- PercentageFeatureSet(CC4, pattern = "^MT-")

CAD1[["percent.mt"]] <- PercentageFeatureSet(CAD1, pattern = "^MT-")
CAD2[["percent.mt"]] <- PercentageFeatureSet(CAD2, pattern = "^MT-")
CAD3[["percent.mt"]] <- PercentageFeatureSet(CAD3, pattern = "^MT-")
CAD4[["percent.mt"]] <- PercentageFeatureSet(CAD4, pattern = "^MT-")

#merge
AD_combined <- merge(HC1,y=list(HC2,HC3,HC4,HC5,HC6,HC7,HC8,
                                HAD1,HAD2,HAD3,HAD4,HAD5,HAD6,HAD7,HAD8,HAD9,
                                CC1,CC2,CC3,CC4,
                                CAD1,CAD2,CAD3,CAD4),
                     add.cell.ids = c("d1","d2","d3","d5","d6"),merge.data = TRUE)

#####################
AD_list <- SplitObject(AD_combined, split.by = "orig.ident")

AD_list <- lapply(X = AD_list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = AD_list)

AD.anchors <- FindIntegrationAnchors(object.list = AD_list, anchor.features = features)

Sys.setenv('R_MAX_VSIZE'=32000000000)
memory.limit(24000)

# this command creates an 'integrated' data assay
AD.combined <- IntegrateData(anchorset = AD.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(AD.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
AD.combined <- ScaleData(AD.combined, verbose = FALSE)
AD.combined <- RunPCA(AD.combined, npcs = 30, verbose = FALSE)
AD.combined <- RunUMAP(AD.combined, reduction = "pca", dims = 1:16)
AD.combined <- FindNeighbors(AD.combined, reduction = "pca", dims = 1:16)
# AD.combined <- FindClusters(AD.combined, resolution = 0.5)
AD.combined <- FindClusters(AD.combined, resolution = 0.15)


# pdf(file = "scomb.umap.pdf",width = 12,height = 8)
# Visualization
DimPlot(AD.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(AD.combined, reduction = "umap")
p1 + p2
# dev.off()

DimPlot(AD.combined, reduction = "umap", split.by = "orig.ident")

DimPlot(AD.combined, reduction = "pca", split.by = "orig.ident")