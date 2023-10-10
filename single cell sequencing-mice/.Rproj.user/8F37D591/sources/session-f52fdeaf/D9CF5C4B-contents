library(Seurat)
library(tidyverse)
library(reticulate)
library(patchwork)
reticulate::py_install(packages = 'umap-learn')

# S0 <- Read10X(data.dir="./S0/filtered_feature_bc_matrix/")
# S0 <- CreateSeuratObject(S0,project = 'M1M3s0',min.cells = 3, min.features = 200)
# S1 <- Read10X(data.dir="./S1/filtered_feature_bc_matrix/")
# S1 <- CreateSeuratObject(S1,project = 'M1M3s1',min.cells = 3, min.features = 200)
# S3 <- Read10X(data.dir="./S3/filtered_feature_bc_matrix/")
# S3 <- CreateSeuratObject(S3,project = 'M1M3s3',min.cells = 3, min.features = 200)
# S7 <- Read10X(data.dir="./S7/filtered_feature_bc_matrix/")
# S7 <- CreateSeuratObject(S7,project = 'M1M3s7',min.cells = 3, min.features = 200)
# S14 <- Read10X(data.dir="./S14/filtered_feature_bc_matrix/")
# S14 <- CreateSeuratObject(S14,project = 'M1M3s14',min.cells = 3, min.features = 200)

# saveRDS('S0.rds')
# saveRDS('S1.rds')
# saveRDS('S3.rds')
# saveRDS('S7.rds')
# saveRDS('S14.rds')

# 1. load data
S0 <- readRDS('S0.rds')
S1 <- readRDS('S1.rds')
S3 <- readRDS('S3.rds')
S7 <- readRDS('S7.rds')
S14 <- readRDS('S14.rds')



S0[['percent.mt']] <- PercentageFeatureSet(S0,pattern = '^mt-')
S1[['percent.mt']] <- PercentageFeatureSet(S1,pattern = '^mt-')
S3[['percent.mt']] <- PercentageFeatureSet(S3,pattern = '^mt-')
S7[['percent.mt']] <- PercentageFeatureSet(S7,pattern = '^mt-')
S14[['percent.mt']] <- PercentageFeatureSet(S14,pattern = '^mt-')


# 2. merge
S_combined <- merge(S0,y = list(S1,S3,S7,S14),
                    add.cell.ids = c('s0','s1','s3','s7','s14'),
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

# save workspace for all integrated samples (0,1,3,7,14)
save.image('S_all.RData')



#  QC plot

pdf(file="integration_plot/QC_plot.pdf",width=18,height=9)
VlnPlot(S.combined,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1 <- FeatureScatter(S.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
# Use "+" from patchwork to replace CombinePlots function from Seurat
plot1 + plot2
dev.off()


# Run the standard workflow for visualization and clustering
S.combined <- ScaleData(S.combined, verbose = FALSE)
S.combined <- RunPCA(S.combined, npcs = 30, verbose = FALSE)
S.combined <- RunUMAP(S.combined, reduction = "pca", dims = 1:10)
S.combined <- FindNeighbors(S.combined, dims = 1:10,k.param = 30,n.trees = 100)
S.combined <- FindClusters(S.combined, resolution = 0.15,n.iter = 10000,n.start = 100)



pdf(file = "integration_plot/scomb.umap.pdf",width = 12,height = 8)
# Visualization
p1 <- DimPlot(S.combined, reduction = "umap", group.by = "orig.ident",label = TRUE)
p2 <- DimPlot(S.combined, reduction = "umap",label = TRUE)
p1 + p2
dev.off()





###### umap plot expand grid

for(dim in 10:20)
{
  S.combined <- RunUMAP(S.combined,reduction = 'pca',dims = 1:dim)
  for(k in 10:30)
  {
    for(reso in seq(0.2,0.4,0.02))
    {
      S.combined <- FindNeighbors(S.combined,
                                  dims = 1:dim,
                                  k.param = k,
                                  n.tree = 100)
      S.combined <- FindClusters(S.combined,
                                 resolution = reso,n.iter = 10000,n.start = 100)
      png(file = paste("integration_plot/expand_grid/scomb_dim",dim,'_k',k,'reso',reso,'_umap.png',sep = ''))
      print(DimPlot(S.combined,reduction = 'umap',label = TRUE))
      dev.off()
    }
  }
}

sink('./integration_plot/expand_grid/platform_info.txt')
version
sessionInfo()

sink()

