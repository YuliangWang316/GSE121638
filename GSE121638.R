library(patchwork)
library(dplyr)
library(Seurat)

T1.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/T1/")
T2.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/T2/")
T3.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/T3/")
P1.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/P1/")
P2.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/P2/")
P3.data<-Read10X("c:/Users/xjmik/Downloads/GSE121638_RAW/singlecell/P3/")

T1.data <- as.data.frame(T1.data)
T2.data <- as.data.frame(T2.data)
T3.data <- as.data.frame(T3.data)
P1.data <- as.data.frame(P1.data)
P2.data <- as.data.frame(P2.data)
P3.data <- as.data.frame(P3.data)

for (i in 1:length(colnames(T1.data))) {
  colnames(T1.data)[i] <- paste(colnames(T1.data)[i],"T1",i,sep = "-")  
}

for (i in 1:length(colnames(T2.data))) {
  colnames(T2.data)[i] <- paste(colnames(T2.data)[i],"T2",i,sep = "-")  
}

for (i in 1:length(colnames(T3.data))) {
  colnames(T3.data)[i] <- paste(colnames(T3.data)[i],"T3",i,sep = "-")  
}

for (i in 1:length(colnames(P1.data))) {
  colnames(P1.data)[i] <- paste(colnames(P1.data)[i],"P1",i,sep = "-")  
}

for (i in 1:length(colnames(P2.data))) {
  colnames(P2.data)[i] <- paste(colnames(P2.data)[i],"P2",i,sep = "-")  
}

for (i in 1:length(colnames(P3.data))) {
  colnames(P3.data)[i] <- paste(colnames(P3.data)[i],"P3",i,sep = "-")  
}


T1.metadata<-data.frame(colnames(T1.data),rep("T1",length(colnames(T1.data))),rep("T",length(colnames(T1.data))))
T2.metadata<-data.frame(colnames(T2.data),rep("T2",length(colnames(T2.data))),rep("T",length(colnames(T2.data))))                                              
T3.metadata<-data.frame(colnames(T3.data),rep("T3",length(colnames(T3.data))),rep("T",length(colnames(T3.data))))
P1.metadata<-data.frame(colnames(P1.data),rep("P1",length(colnames(P1.data))),rep("P",length(colnames(P1.data))))
P2.metadata<-data.frame(colnames(P2.data),rep("P2",length(colnames(P2.data))),rep("P",length(colnames(P2.data))))
P3.metadata<-data.frame(colnames(P3.data),rep("P3",length(colnames(P3.data))),rep("P",length(colnames(P3.data))))

colnames(T1.metadata)<-c("barcode","group1","group2")
colnames(T2.metadata)<-c("barcode","group1","group2")
colnames(T3.metadata)<-c("barcode","group1","group2")
colnames(P1.metadata)<-c("barcode","group1","group2")
colnames(P2.metadata)<-c("barcode","group1","group2")
colnames(P3.metadata)<-c("barcode","group1","group2")

pbmc.metadata<-rbind(T1.metadata,T2.metadata,T3.metadata,P1.metadata,P2.metadata,P3.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(T1.data,T2.data,T3.data,P1.data,P2.data,P3.data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200,meta.data = pbmc.metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc)
FeaturePlot(pbmc,features = "FOXP3")
VlnPlot(pbmc,features = "FOXP3",sort = TRUE,pt.size = 0)
Treg<-subset(pbmc,idents = "12")
Idents(Treg)<-Treg@meta.data$group2
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggpubr)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
FeaturePlot(Treg,features = "JMJD1C",split.by = "group2",order = TRUE)

rownames(T1.metadata)<-T1.metadata[,1]
rownames(T2.metadata)<-T2.metadata[,1]
rownames(T3.metadata)<-T3.metadata[,1]
rownames(P1.metadata)<-P1.metadata[,1]
rownames(P2.metadata)<-P2.metadata[,1]
rownames(P3.metadata)<-P3.metadata[,1]

T1 <- CreateSeuratObject(counts = T1.data, project = "IMMUNE_T1", min.cells = 5,meta.data = T1.metadata)
T1$type <- "T1"
T1 <- subset(T1, subset = nFeature_RNA > 500)
T1 <- NormalizeData(T1, verbose = FALSE)
T1 <- FindVariableFeatures(T1, selection.method = "vst", nfeatures = 2000)

T2 <- CreateSeuratObject(counts = T2.data, project = "IMMUNE_T2", min.cells = 5,meta.data = T2.metadata)
T2$type <- "T2"
T2 <- subset(T2, subset = nFeature_RNA > 500)
T2 <- NormalizeData(T2, verbose = FALSE)
T2 <- FindVariableFeatures(T2, selection.method = "vst", nfeatures = 2000)

T3 <- CreateSeuratObject(counts = T3.data, project = "IMMUNE_T3", min.cells = 5,meta.data = T3.metadata)
T3$type <- "T3"
T3 <- subset(T3, subset = nFeature_RNA > 500)
T3 <- NormalizeData(T3, verbose = FALSE)
T3 <- FindVariableFeatures(T3, selection.method = "vst", nfeatures = 2000)

P1 <- CreateSeuratObject(counts = P1.data, project = "IMMUNE_P1", min.cells = 5,meta.data = P1.metadata)
P1$type <- "P1"
P1 <- subset(P1, subset = nFeature_RNA > 500)
P1 <- NormalizeData(P1, verbose = FALSE)
P1 <- FindVariableFeatures(P1, selection.method = "vst", nfeatures = 2000)

P2 <- CreateSeuratObject(counts = P2.data, project = "IMMUNE_P2", min.cells = 5,meta.data = P2.metadata)
P2$type <- "P2"
P2 <- subset(P2, subset = nFeature_RNA > 500)
P2 <- NormalizeData(P2, verbose = FALSE)
P2 <- FindVariableFeatures(P2, selection.method = "vst", nfeatures = 2000)

P3 <- CreateSeuratObject(counts = P3.data, project = "IMMUNE_P3", min.cells = 5,meta.data = P3.metadata)
P3$type <- "P3"
P3 <- subset(P3, subset = nFeature_RNA > 500)
P3 <- NormalizeData(P3, verbose = FALSE)
P3 <- FindVariableFeatures(P3, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(T1, T2,T3,P1,P2,P3), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")
FeaturePlot(immune.combined,features = "FOXP3")
VlnPlot(immune.combined,features = "FOXP3",sort = TRUE,pt.size = 0,assay = "RNA")
Treg<-subset(immune.combined,idents = "11")
Idents(Treg)<-Treg@meta.data$group2
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggpubr)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE,pt.size = 0)+stat_compare_means()
FeaturePlot(Treg,features = "JMJD1C",split.by = "group2",order = TRUE)
FeaturePlot(Treg,features = "JMJD1C",split.by = "group2")

JMJD1C<-FetchData(Treg,vars = "JMJD1C")
JMJD1C_metadata<-data.frame(rownames(Treg@meta.data),Treg@meta.data$group2)
JMJD1C_new<-cbind(JMJD1C,JMJD1C_metadata[,2])
write.table(JMJD1C_new,file = "GSE121638.txt",sep = "\t")
