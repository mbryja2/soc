
library(dplyr)
library(Seurat)
library(patchwork)

rds<- file.choose()
tmcrds <- readRDS(rds)

markers1<-FindMarkers(tmcrds, ident.1 = "Basal 1", min.pct = 0.25)


DimPlot(tmcrds, dims = 1:2, reduction = "umap")


VlnPlot(tmcrds,features = "TXNIP")
head(markers1)


tmcrds.markers<-FindAllMarkers(tmcrds,only.pos = TRUE, min.pct = 0.25)

topmarkers <-tmcrds.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )
