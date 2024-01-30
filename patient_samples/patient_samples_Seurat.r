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

sample_OSCC_17.data<-Read10X(data.dir = "OSCC_17/outs/filtered_feature_bc_matrix")
sample_OSCC_17 = CreateSeuratObject(counts = sample_OSCC_17.data, project = "Sample_OSCC_17", min.cells = 3, min.features = 200)
sample_OSCC_17[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_17, pattern = "^MT-")
sample_OSCC_12.data<-Read10X(data.dir = "OSCC_12/outs/filtered_feature_bc_matrix")
sample_OSCC_12 = CreateSeuratObject(counts = sample_OSCC_12.data, project = "Sample_OSCC_12", min.cells = 3, min.features = 200)
sample_OSCC_12[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_12, pattern = "^MT-")
sample_OSCC_13.data<-Read10X(data.dir = "OSCC_13/outs/filtered_feature_bc_matrix")
sample_OSCC_13 = CreateSeuratObject(counts = sample_OSCC_13.data, project = "Sample_OSCC_13", min.cells = 3, min.features = 200)
sample_OSCC_13[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_13, pattern = "^MT-")
sample_OSCC_14.data<-Read10X(data.dir = "OSCC_14/outs/filtered_feature_bc_matrix")
sample_OSCC_14 = CreateSeuratObject(counts = sample_OSCC_14.data, project = "Sample_OSCC_14", min.cells = 3, min.features = 200)
sample_OSCC_14[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_14, pattern = "^MT-")
sample_OSCC_15.data<-Read10X(data.dir = "OSCC_15/outs/filtered_feature_bc_matrix")
sample_OSCC_15 = CreateSeuratObject(counts = sample_OSCC_15.data, project = "Sample_OSCC_15", min.cells = 3, min.features = 200)
sample_OSCC_15[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_15, pattern = "^MT-")
sample_OSCC_11.data<-Read10X(data.dir = "OSCC_11/outs/filtered_feature_bc_matrix")
sample_OSCC_11 = CreateSeuratObject(counts = sample_OSCC_11.data, project = "Sample_OSCC_11", min.cells = 3, min.features = 200)
sample_OSCC_11[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_11, pattern = "^MT-")
sample_OSCC_16.data<-Read10X(data.dir = "OSCC_16/outs/filtered_feature_bc_matrix")
sample_OSCC_16 = CreateSeuratObject(counts = sample_OSCC_16.data, project = "Sample_OSCC_16", min.cells = 3, min.features = 200)
sample_OSCC_16[["percent.mt"]] <- PercentageFeatureSet(sample_OSCC_16, pattern = "^MT-")
# merge, cluster and Harmony integration
sample.headneck <- merge(sample_OSCC_17, y = c(sample_OSCC_12,sample_OSCC_13,sample_OSCC_14,sample_OSCC_15,sample_OSCC_11,sample_OSCC_16), add.cell.ids = c('Sample_OSCC_17','Sample_OSCC_12','Sample_OSCC_13_T','Sample_OSCC_14_T','Sample_OSCC_15_T','Sample_OSCC_11',"Sample_OSCC_16"), project = "SAMPLE.INTEGRATED")

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
