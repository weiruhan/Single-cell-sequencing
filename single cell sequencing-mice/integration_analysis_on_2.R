# Integrate S0 and S1
library(Seurat)
library(tidyverse)
library(reticulate)
library(patchwork)
reticulate::py_install(packages = 'umap-learn')

# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap')
# remotes::install_github("metaOmics/MetaDE")

library(BiocManager)
library(multtest)
library(metap)
# library(metaDE)

S0 <- readRDS('S0.rds')
S1 <- readRDS('S1.rds')

S0[['percent.mt']] <- PercentageFeatureSet(S0,pattern = '^mt-')
S1[['percent.mt']] <- PercentageFeatureSet(S1,pattern = '^mt-')

# 2. merge
S_combined <- merge(S0,y = S1,
                    add.cell.ids = c('s0','s1'),
                    merge.data = TRUE)
S_list <- SplitObject(S_combined,split.by = 'orig.ident')

S_list <- lapply(X = S_list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 3000 & percent.mt < 1.5)
  x <- NormalizeData(x, normalization.method = "LogNormalize",
                     scale.factor = 10000)
  # x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = S_list)

S.anchors <- FindIntegrationAnchors(object.list = S_list, 
                                    anchor.features = features)


# this command creates an 'integrated' data assay
S.combined <- IntegrateData(anchorset = S.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(S.combined) <- "integrated"




# Run the standard workflow for visualization and clustering
S.combined <- ScaleData(S.combined, verbose = FALSE)
S.combined <- RunPCA(S.combined, npcs = 30, verbose = FALSE)
S.combined <- RunUMAP(S.combined, reduction = "pca", dims = 1:30)
S.combined <- FindNeighbors(S.combined, dims = 1:30,k.param = 30,n.trees = 100)
S.combined <- FindClusters(S.combined, resolution = 0.3,n.iter = 10000,n.start = 100)

# save.image('S01.RData')


pdf(file = "integration_plot/scomb_S0S1.umap.pdf",width = 12,height = 8)
# Visualization
p1 <- DimPlot(S.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(S.combined, reduction = "umap",label = TRUE)
p3 <- DimPlot(S.combined, reduction = 'umap',split.by = 'orig.ident')
p1 + p2
p3
dev.off()

# Integrate S0 and S14

S0 <- readRDS('S0.rds')
S14 <- readRDS('S14.rds')

S0[['percent.mt']] <- PercentageFeatureSet(S0,pattern = '^mt-')
S14[['percent.mt']] <- PercentageFeatureSet(S14,pattern = '^mt-')

# 2. merge
S_combined <- merge(S0,y = S14,
                    add.cell.ids = c('s0','s14'),
                    merge.data = TRUE)
S_list <- SplitObject(S_combined,split.by = 'orig.ident')

S_list <- lapply(X = S_list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 3000 & percent.mt < 1.5)
  x <- NormalizeData(x, normalization.method = "LogNormalize",
                     scale.factor = 10000)
  # x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = S_list)

S.anchors <- FindIntegrationAnchors(object.list = S_list, 
                                    anchor.features = features)


# this command creates an 'integrated' data assay
S.combined <- IntegrateData(anchorset = S.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(S.combined) <- "integrated"




# Run the standard workflow for visualization and clustering
S.combined <- ScaleData(S.combined, verbose = FALSE)
S.combined <- RunPCA(S.combined, npcs = 30, verbose = FALSE)
S.combined <- RunUMAP(S.combined, reduction = "pca", dims = 1:30)
S.combined <- FindNeighbors(S.combined, dims = 1:30,k.param = 30,n.trees = 100)
S.combined <- FindClusters(S.combined, resolution = 0.3,n.iter = 10000,n.start = 100)

# save.image('S014.RData')

pdf(file = "integration_plot/scomb_S0S14.umap.pdf",width = 12,height = 8)
# Visualization
p1 <- DimPlot(S.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(S.combined, reduction = "umap",label = TRUE)
p3 <- DimPlot(S.combined, reduction = 'umap',split.by = 'orig.ident')
p1 + p2
p3
dev.off()