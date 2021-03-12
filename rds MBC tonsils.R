
library(dplyr)
library(Seurat)
library(patchwork)

rds<- file.choose()
tmcrds <- readRDS(rds)

markersbasal1x<-FindMarkers(tmcrds, ident.1 = "Basal 1", min.pct = 0.25)
markersfcrl<- FindMarkers(tmcrds, ident.1 = "FCRL2/3high MBC" )
Idents(tmcrds)<-tmcrds$Donor

?FindMarkers

memory.limit(size = 200000)

DimPlot(tmcrds, dims = 1:2, reduction = "umap", label=TRUE)

        
VlnPlot(tmcrds,features = "RPS27")
head(markers1)


tmcrds.markers<-FindAllMarkers(tmcrds,only.pos = TRUE, min.pct = 0.25)

topmarkers <-tmcrds.markers%>% group_by(cluster)%>% top_n(n=2, wt= avg_logFC )

markers <- tmcrds.markers%>%group_by(cluster)

markersbasal1 <- subset(markers, cluster %in% "Basal 1")
markersactivated1 <- subset(markers, cluster %in% "Activated 1")
markersfcrl <- subset(markers, cluster %in% "FCRL2/3high MBC")



