#!/usr/bin/env Rscript
library(harmony)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(msigdbr)
library(cowplot)
library(dplyr)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
library("enrichplot")
library(ggupset)
library(gridExtra)
library(pheatmap)
options(bitmapType = 'cairo')
knitr::opts_chunk$set(dev="CairoPNG")
library(org.Hs.eg.db)

sample_9347.data<-Read10X(data.dir = "9347/outs/filtered_feature_bc_matrix")
sample_9347 = CreateSeuratObject(counts = sample_9347.data, project = "Sample_9347", min.cells = 3, min.features = 200)
sample_9347[["percent.mt"]] <- PercentageFeatureSet(sample_9347, pattern = "^MT-")
sample_9398.data<-Read10X(data.dir = "9398/outs/filtered_feature_bc_matrix")
sample_9398 = CreateSeuratObject(counts = sample_9398.data, project = "Sample_9398", min.cells = 3, min.features = 200)
sample_9398[["percent.mt"]] <- PercentageFeatureSet(sample_9398, pattern = "^MT-")
sample_9236.data<-Read10X(data.dir = "9236/outs/filtered_feature_bc_matrix")
sample_9236 = CreateSeuratObject(counts = sample_9236.data, project = "Sample_9236", min.cells = 3, min.features = 200)
sample_9236[["percent.mt"]] <- PercentageFeatureSet(sample_9236, pattern = "^MT-")
sample_9237.data<-Read10X(data.dir = "9237/outs/filtered_feature_bc_matrix")
sample_9237 = CreateSeuratObject(counts = sample_9237.data, project = "Sample_9237", min.cells = 3, min.features = 200)
sample_9237[["percent.mt"]] <- PercentageFeatureSet(sample_9237, pattern = "^MT-")
sample_9218.data<-Read10X(data.dir = "9218/outs/filtered_feature_bc_matrix")
sample_9218 = CreateSeuratObject(counts = sample_9218.data, project = "Sample_9218", min.cells = 3, min.features = 200)
sample_9218[["percent.mt"]] <- PercentageFeatureSet(sample_9218, pattern = "^MT-")
sample_BM320030.data<-Read10X(data.dir = "BM320030/outs/filtered_feature_bc_matrix")
sample_BM320030 = CreateSeuratObject(counts = sample_BM320030.data, project = "Sample_BM320030", min.cells = 3, min.features = 200)
sample_BM320030[["percent.mt"]] <- PercentageFeatureSet(sample_BM320030, pattern = "^MT-")
sample_BM319435.data<-Read10X(data.dir = "BM319435/outs/filtered_feature_bc_matrix")
sample_BM319435 = CreateSeuratObject(counts = sample_BM319435.data, project = "Sample_BM319435", min.cells = 3, min.features = 200)
sample_BM319435[["percent.mt"]] <- PercentageFeatureSet(sample_BM319435, pattern = "^MT-")
# merge, cluster and Harmony integration
sample.headneck <- merge(sample_9347, y = c(sample_9398,sample_9236,sample_9237,sample_9218,sample_BM320030,sample_BM319435), add.cell.ids = c('Sample_9347','Sample_9398','Sample_9236','Sample_9237','Sample_9218','Sample_BM320030',"Sample_BM319435"), project = "SAMPLE.INTEGRATED")

VlnPlot(sample.headneck, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sample.headneck <- NormalizeData(object = sample.headneck, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssay(sample.headneck,assay = "RNA")
sample.headneck <- FindVariableFeatures(object = sample.headneck, selection.method = "vst", nfeatures = 5000)
top20 <- head(x = VariableFeatures(object = sample.headneck), 20)
plot1 <- VariableFeaturePlot(object = sample.headneck)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1+plot2
all.genes <- rownames(sample.headneck)
sample.headneck<- ScaleData(object = sample.headneck,features = all.genes)
sample.headneck <- RunPCA(object = sample.headneck,pc.genes = VariableFeatures(sample.headneck))
ElbowPlot(sample.headneck)
sample.headneck <- RunHarmony(sample.headneck, group.by.vars="orig.ident",assay.use='RNA')
names(sample.headneck@reductions)
sample.headneck <- RunUMAP(sample.headneck,  dims = 1:20, 
                     reduction = "harmony",seed.use=111)
DimPlot(sample.headneck,reduction = "umap",label=T ) 
DimPlot(sample.headneck,reduction = "umap",label=F ) + ggtitle("Integrated")+theme(plot.title = element_text(hjust = 0.5))

sample.headneck <- FindNeighbors(sample.headneck, reduction = "harmony",dims = 1:20)
sample.headneck <- FindClusters(sample.headneck, resolution = 0.5)
table(sample.headneck@meta.data$seurat_clusters)
DimPlot(sample.headneck,reduction = "umap",label=T)  
DimPlot(sample.headneck,reduction = "umap",label=T,
        group.by = 'orig.ident') 

# SingleR annotation
ref <- HumanPrimaryCellAtlasData()
seuratObj_annot <- as.SingleCellExperiment(sample.headneck)
library(SingleR)
pred <- SingleR(test=seuratObj_annot, ref=ref, labels=ref$label.fine)
head(pred)
plotScoreHeatmap(pred)
tab <- table(Assigned=pred$pruned.labels, Cluster=seuratObj_annot@colData$seurat_clusters)
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.

pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
pred2 <- SingleR(test=seuratObj_annot, ref=ref, cluster=seuratObj_annot@colData$seurat_clusters, labels=ref$label.fine)
sample.headneck.backup = sample.headneck
sample.headneck@meta.data$cell.type.fine = sample.headneck@meta.data$seurat_clusters
sample.headneck[["SingleR.cluster.labels"]] <- 
        pred2$labels[match(sample.headneck[[]][["seurat_clusters"]], rownames(pred2))]

Idents(sample.headneck) <- "SingleR.cluster.labels"
sample.headneck <- RenameIdents(sample.headneck, 
    'Epithelial_cells:bronchial' = "Epithelial_cells",
    'Epithelial_cells:bladder' = "Epithelial_cells"
    )

# CopyKAT prediction
sample.headneck.exp.rawdata <- as.matrix(sample.headneck@assays$RNA@counts)
sample.headneck.copykat.test <- copykat(rawmat=sample.headneck.exp.rawdata, id.type="S", ngene.chr=3, win.size=25, KS.cut=0.1, sam.name="sample.headneck", distance="euclidean", norm.cell.names="", n.cores=34,output.seg="FLASE")
# Add CopyKAT metadata
sample.headneck_copykat_prediction_csv = 'sample.headneck_copykat_prediction.txt'
sample.headneck_copykat_prediction<-read.csv(sample.headneck_copykat_prediction_csv,sep='\t',header=TRUE,row.names=1)
sample.headneck<-AddMetaData(sample.headneck, sample.headneck_copykat_prediction)

# add pathogen UMI metadata
umi_table_csv = 'csv_novami_mix_dedup_rename.csv'
umi_table<-read.csv(umi_table_csv,sep=',',header=TRUE,row.names = 1)
umi_table[is.na(umi_table)] <- 0
umi_table$Total <- rowSums(umi_table)
umi_table[umi_table==0] <- NA
sample.headneck<-AddMetaData(sample.headneck, umi_table)
#saveRDS(sample.headneck, file = sample.headneck.rds)

