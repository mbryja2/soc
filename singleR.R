library(dplyr)
library(Seurat)
library(patchwork)
library(BiocManager)
library(SingleR)
library(SingleCellExperiment)
library(viridis)
library(pheatmap)
library(installr)

install.packages("SingleR")
install.packages("Biocmanager")

memory.limit(size = 200000 )

BiocManager::install("SingleCellExperiment")

tmcse <- as.SingleCellExperiment(tmc)
tmcrdsse <-as.SingleCellExperiment (tmcrds)

cllU<-as.SingleCellExperiment(cllusT)
cllM<-as.SingleCellExperiment(cllmsT)



tmcn<- GetAssayData(tmc)
tmcrdsn <- GetAssayData(tmcrds)

#com<- SingleR(test = tmcse, ref = tmcrdsse,assay.type.test = "logcounts" , assay.type.ref = "logcounts", labels = tmcrdsse$CellType)
comc<- SingleR(test = tmcse,clusters = tmc$seurat_clusters, ref = tmcrdsse,assay.type.test = "logcounts" , assay.type.ref = "logcounts", labels = tmcrdsse$CellType)

prdel<-SingleR(test = tmcse, ref = tmcrdsse,assay.type.test = "logcounts" , assay.type.ref = "logcounts", labels = tmcrdsse$CellType)
plotScoreHeatmap(prdel)

?SingleR

table(com$pruned.labels)

plotScoreHeatmap(comc)

plotDeltaDistribution(com)


topn<-top_n(com$scores["CD21low2"])
     
anncllM<- SingleR(test = cllM, ref = tmcrdsse, assay.type.test = "logcounts" , assay.type.ref = "logcounts", labels = tmcrdsse$CellType)     

?DimPlot
plotScoreHeatmap(anncllM)

anncllM
     
anncllU<- SingleR(test = cllU, ref = tmcrdsse, assay.type.test = "logcounts" , assay.type.ref = "logcounts", labels = tmcrdsse$CellType)     

plotScoreHeatmap(anncllU)

com

tmc[["SingleR.labels"]]<- com$labels

cllusT[["SingleR.labels"]]<- anncllU$labels

cllmsT[["SingleR.labels"]]<- anncllM$labels

