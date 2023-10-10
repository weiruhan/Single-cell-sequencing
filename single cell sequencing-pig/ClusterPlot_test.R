# Use 15 PCs to run TSNE
# pbmc<-RunTSNE(pbmc,dims=1:15)
# pbmc <- RunUMAP(pbmc,dims = 1:15)
# Use KNN to run do clustering
pbmc <- FindNeighbors(pbmc, dims = 1:15,k.param = 20,n.trees = 100)
pbmc <- FindClusters(pbmc, resolution = 0.4,n.iter = 10000,start = 100)

DimPlot(pbmc,reduction = 'umap',label = TRUE)

for(um in 10:15)
{
  pbmc <- RunUMAP(pbmc,dims = 1:um)
  for(k in c(10,20,30,40,50))
  {
    for(reso in seq(0.2,0.4,0.02))
    {
      plot_name <- paste('plot/umap_',um,'pc_',k,'kpar_',reso,'reso','.jpeg',sep = '')
      pbmc <- FindNeighbors(pbmc, dims = 1:um,k.param = k,n.trees = 100)
      pbmc <- FindClusters(pbmc, resolution = reso,n.iter = 10000,start = 100)
      figure <- DimPlot(pbmc,reduction = 'umap',label = TRUE)
      ggsave(figure,file = plot_name)
    }
  }
}