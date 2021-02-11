library(dplyr)
library(Seurat)
library(patchwork)

#sorted memory cells, tonsils
data_dir <- "D:/Users/Jipro/Downloads/BCP004_MBC_5GEX_filtered_feature_bc_matrix/"
#list.files(data_dir)

tmc_data <- Read10X("D:/Users/Jipro/Downloads/BCP004_MBC_5GEX_filtered_feature_bc_matrix/")


tmc <-CreateSeuratObject(counts = tmc_data )
tmc[["percent.mt"]]<- PercentageFeatureSet(tmc, pattern = "^MT-")

    VlnPlot(tmc, features = c("percent.mt","nCount_RNA", "nFeature_RNA"), ncol = 3)
    VlnPlot(tmc, features = "percent.mt")
    
plot1 <- FeatureScatter(tmc,"nCount_RNA", "nFeature_RNA")
plot2 <- FeatureScatter(tmc,"nCount_RNA", "percent.mt")
#plot1
#plot2

tmc <- subset(tmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2800 & percent.mt < 7)

tmc<- NormalizeData(tmc)

tmc<- FindVariableFeatures(tmc,selection.method = "vst" )
top10 <- head(VariableFeatures(tmc), 20)

plot1v <-VariableFeaturePlot(tmc)
plot2v <-LabelPoints(plot1v, points = top10)
#plot2v

all.genes <- rownames(tmc)
tmc <- ScaleData(tmc, features = all.genes)

tmc<- RunPCA(tmc, features = VariableFeatures(tmc))

#print(tmc[["pca"]], dims = 1)

VizDimLoadings(tmc, dims = 1)

DimPlot(tmc, dims = c(2,6))

DimHeatmap(tmc,dims = 1,cells = 500)

tmc <- JackStraw(tmc, num.replicate = 100)
tmc <- ScoreJackStraw(tmc, dims = 1:20)

#JackStrawPlot(tmc, dims = 1:20)


ElbowPlot(tmc)  

tmc <- FindNeighbors(tmc, dims = 1:10) 
tmc<- FindClusters(tmc, resolution = 1) 
 
head(Idents(tmc))

tmc<- RunUMAP(tmc, dims = 1:10)

DimPlot(tmc, dims =1:2, reduction = "umap")

cluster1.markers <- FindMarkers(tmc, ident.1 = 1, min.pct = 0.25) 
head(cluster1.markers, n=30) 
cluster2.markers <- FindMarkers(tmc, ident.1 = 2, min.pct = 0.25) 


cluster3.markers <- FindMarkers(tmc, ident.1= 1, ident.2 = 2, min.pct = 0.25)
#head(cluster3.markers)


tmc.markers<-FindAllMarkers(tmc,only.pos = TRUE, min.pct = 0.25)

tmc.topmarkers <-tmc.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )

#VlnPlot(tmc, features = c("ZFP36L2", "JUNB", "FCER2"))

#FeaturePlot(tmc, features = c("ZFP36L2", "JUNB"))





top10m <- tmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

head(cluster1.markers)
                                                      




