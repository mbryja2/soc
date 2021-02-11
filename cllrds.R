library(dplyr)
library(Seurat)
library(patchwork)

install.packages('BiocManager')
BiocManager::install('limma')

rds<- file.choose()


cllrds<- readRDS(rds)

#cllrds[["percent.mt"]]<- PercentageFeatureSet(cllrds, pattern = "^MT-")

#cllrds <- subset(cllrds,subset = nFeature_RNA > 200 & nFeature_RNA < 2800 & percent.mt < 20)

#cllrds<- NormalizeData(cllrds)

#cllrds<- ScaleData(cllrds, features = row.names(cllrds))

#cllrds <- FindVariableFeatures(cllrds, selection.method = "vst")

#cllrds <- RunPCA(cllrds, features = VariableFeatures(cllrds))

#cllrds<- RunUMAP(cllrds, dims = 1:5)
cllrds<- FindNeighbors(cllrds)
cllrds <- FindClusters(cllrds,resolution = 0.1)



DimPlot(cllrds, reduction = "umap" )

cllrds.markers<- FindAllMarkers(cllrds, min.pct = 0.25)

topcllmrk<-cllrds.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )

