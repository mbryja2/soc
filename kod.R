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

?VlnPlot
    
plot1 <- FeatureScatter(tmc,"nCount_RNA", "nFeature_RNA")
plot2 <- FeatureScatter(tmc,"nCount_RNA", "percent.mt")
#plot1
#plot2

tmc <- subset(tmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2800 & percent.mt < 7)

tmc<- NormalizeData(tmc)

tmc<- FindVariableFeatures(tmc,selection.method = "vst" )
top10 <- head(VariableFeatures(tmc), 5)

plot1v <-VariableFeaturePlot(tmc)
plot2v <-LabelPoints(plot1v, points = top10)
plot2v

#all Ig V, D, J genes (extracted using the regular expression [regex] “IG[HKL][VDJ]”),
VariableFeatures(tmc) <- VariableFeatures(tmc)[-grep("IG[HKL][VDJ]", VariableFeatures(tmc))]
#Ig constant genes (IGHM, IGHD, IGHE, IGHA[1-2], IGHG[1-4], IGKC, IGLC[1-7], and AC233755.1 [which encodes IGHV4-38-2]),
VariableFeatures(tmc) <- VariableFeatures(tmc)[-grep("IGHA[1-2]", VariableFeatures(tmc))]
VariableFeatures(tmc) <- VariableFeatures(tmc)[-grep("IGHG[1-4]", VariableFeatures(tmc))]
VariableFeatures(tmc) <- VariableFeatures(tmc)[-grep("IGLC[1-7]", VariableFeatures(tmc))]

VariableFeatures(tmc) <- VariableFeatures(tmc)[VariableFeatures(tmc)!="IGHM"]
VariableFeatures(tmc) <- VariableFeatures(tmc)[VariableFeatures(tmc)!="IGHD"]
VariableFeatures(tmc) <- VariableFeatures(tmc)[VariableFeatures(tmc)!="IGHE"]
VariableFeatures(tmc) <- VariableFeatures(tmc)[VariableFeatures(tmc)!="IGKC"]
VariableFeatures(tmc) <- VariableFeatures(tmc)[VariableFeatures(tmc)!="AC233755.1"]
#T-cell receptor genes (regex “TR[ABGD][CV]”).
VariableFeatures(tmc) <- VariableFeatures(tmc)[-grep("TR[ABGD][CV]", VariableFeatures(tmc))]


all.genes <- rownames(tmc)
tmc <- ScaleData(tmc, features = all.genes)

tmc<- RunPCA(tmc, features = VariableFeatures(tmc))

#print(tmc[["pca"]], dims = 1)

VizDimLoadings(tmc, dims = 1)

DimPlot(tmc, dims = c(2,6))

DimHeatmap(tmc,dims = 1,cells = 500)

#tmc <- JackStraw(tmc, num.replicate = 100)
#tmc <- ScoreJackStraw(tmc, dims = 1:20)

#JackStrawPlot(tmc, dims = 1:20)


ElbowPlot(tmc)  

tmc <- FindNeighbors(tmc, dims = 1:10) 
tmc<- FindClusters(tmc, resolution = 1) 
 
head(Idents(tmc))

tmc<- RunUMAP(tmc, dims = 1:10)

dimplot<- DimPlot(tmc, dims =1:2, reduction = "umap", label = TRUE)

dimplot





cluster1.markers <- FindMarkers(tmc, ident.1 = 1, min.pct = 0.25) 
head(cluster1.markers, n=5) 
cluster2.markers <- FindMarkers(tmc, ident.1 = 2, min.pct = 0.25) 


cluster3.markers <- FindMarkers(tmc, ident.1= 1, ident.2 = 2, min.pct = 0.25)
#head(cluster3.markers)


tmc.markers<-FindAllMarkers(tmc,only.pos = TRUE, min.pct = 0.25)

tmc.topmarkers <-tmc.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )

#VlnPlot(tmc, features = c("ZFP36L2", "JUNB", "FCER2"))

FeaturePlot(tmcrds, features = c("CD5"))

FeaturePlot(tmc, features = c("TSC22D3", "KLF2"))
tmc.topmarkers

DimPlot(tmc, group.by = "SingleR.labels", label = TRUE)

tmc.markers


table(tmc$seurat_clusters, tmc$SingleR.labels)

tmc$
  
top10m <- tmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10m

head(cluster1.markers)
                                                      


cluster0<-subset(tmc, seurat_clusters %in% 0)
cluster1<-subset(tmc, seurat_clusters %in% 1)
cluster2<-subset(tmc, seurat_clusters %in% 2)
cluster3<-subset(tmc, seurat_clusters %in% 3)
cluster4<-subset(tmc, seurat_clusters %in% 4)
cluster5<-subset(tmc, seurat_clusters %in% 5)
cluster6<-subset(tmc, seurat_clusters %in% 6)
cluster7<-subset(tmc, seurat_clusters %in% 7)
cluster8<-subset(tmc, seurat_clusters %in% 8)
cluster9<-subset(tmc, seurat_clusters %in% 9)
cluster10<-subset(tmc, seurat_clusters %in% 10)
cluster11<-subset(tmc, seurat_clusters %in% 11)

table0<- cluster0[[c("seurat_clusters","SingleR.labels")]]
table1<- cluster1[[c("seurat_clusters","SingleR.labels")]]
table2<- cluster2[[c("seurat_clusters","SingleR.labels")]]
table3<- cluster3[[c("seurat_clusters","SingleR.labels")]]
table4<- cluster4[[c("seurat_clusters","SingleR.labels")]]
table5<- cluster5[[c("seurat_clusters","SingleR.labels")]]
table6<- cluster6[[c("seurat_clusters","SingleR.labels")]]
table7<- cluster7[[c("seurat_clusters","SingleR.labels")]]
table8<- cluster8[[c("seurat_clusters","SingleR.labels")]]
table9<- cluster9[[c("seurat_clusters","SingleR.labels")]]
table10<- cluster10[[c("seurat_clusters","SingleR.labels")]]
table11<- cluster11[[c("seurat_clusters","SingleR.labels")]]

head(table10)

nrow(table10$SingleR.labels == "preGC MBC")

