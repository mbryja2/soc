library(dplyr)
library(Seurat)
library(patchwork)
library(BiocManager)
library(SingleR)

install.packages("BiocManager")
BiocManager::install("limma")


rds<- file.choose()


cllrds<- readRDS(rds)

cllrdsU<- subset(cllrds, HTO_classification %in% c("1892-4C-0h", "1892-RT-0h"))
cllrdsM<- subset(cllrds, HTO_classification %in% c("1472-4C-0h", "1472-RT-0h"))

cllmsT<- subset(cllrdsM, seurat_clusters %in% c("1","5"))
cllusT <- subset(cllrdsU, seurat_clusters %in% c("0", "4","3"))                                              

cllrdsU$seurat_clusters

#cllrds[["percent.mt"]]<- PercentageFeatureSet(cllrds, pattern = "^MT-")

#cllrds <- subset(cllrds,subset = nFeature_RNA > 200 & nFeature_RNA < 2800 & percent.mt < 20)

#cllrds<- NormalizeData(cllrds)

#cllrds<- ScaleData(cllrds, features = row.names(cllrds))

#cllrds <- FindVariableFeatures(cllrds, selection.method = "vst")

#cllrds <- RunPCA(cllrds, features = VariableFeatures(cllrds))

#cllrds<- RunUMAP(cllrds, dims = 1:5)
cllrds<- FindNeighbors(cllrds)
cllrds <- FindClusters(cllrds,resolution = 0.1)

cllUBasal1 <- subset(cllrdsU, SingleR.labels %in% c("Basal 1")  )
cllUActivated1 <- subset(cllrdsU, SingleR.labels %in% c("Activated 1")  )
cllUFCRL23highMBC <- subset(cllrdsU, SingleR.labels %in% c("FCRL2/3high MBC")  )
cllUCD21low2 <- subset(cllrdsU, SingleR.labels %in% c("CD21low 2")  )
cllUpreGCMBC <- subset(cllrdsU, SingleR.labels %in% c("preGC MBC")  )
cllUCD21low1 <- subset(cllrdsU, SingleR.labels %in% c("CD21low 1")  )

nubasal1 <- ncol(cllUBasal1)
nuA1<- ncol(cllUActivated1)
nuCD1<- ncol(cllUCD21low1)
nuCD2<- ncol(cllUCD21low2)
nupreGC<- ncol(cllUpreGCMBC)
nuF23<- ncol(cllUFCRL23highMBC)

tableU <-data.frame(row.names =  c("nubasal1", "nuA1", "nuCD1", "nuCD2", "nupreGC", "nuF23"), c(nubasal1, nuA1, nuCD1, nuCD2, nupreGC, nuF23) )


tableU

cllMBasal1 <- subset(cllrdsM, SingleR.labels %in% c("Basal 1")  )
cllMActivated1 <- subset(cllrdsM, SingleR.labels %in% c("Activated 1")  )
cllMFCRL23highMBC <- subset(cllrdsM, SingleR.labels %in% c("FCRL2/3high MBC")  )
cllMCD21low2 <- subset(cllrdsM, SingleR.labels %in% c("CD21low 2")  )
cllMpreGCMBC <- subset(cllrdsM, SingleR.labels %in% c("preGC MBC")  )
cllMCD21low1 <- subset(cllrdsM, SingleR.labels %in% c("CD21low 1")  )

nmbasal1 <- ncol(cllMBasal1)
nmA1<- ncol(cllMActivated1)
nmCD1<- ncol(cllMCD21low1)
nmCD2<- ncol(cllMCD21low2)
nmpreGC<- ncol(cllMpreGCMBC)
nmF23<- ncol(cllMFCRL23highMBC)



tableM <- data.frame(row.names =  c("nmbasal1", "nmA1", "nmCD1", "nmCD2", "nmpreGC", "nmF23"), c(nmbasal1, nmA1, nmCD1, nmCD2, nmpreGC, nmF23) )

tableM

head(cllrdsU$SingleR.labels)

DimPlot(cllmsT, reduction = "umap", label = TRUE )

DimPlot(cllusT, label = TRUE)
DimPlot(cllmsT, label = TRUE)

DimPlot()
cllusT$




cllrds.markers<- FindAllMarkers(cllrds, min.pct = 0.25)

topcllmrk<-cllrds.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )

table(cllrds$HTO_classification)
