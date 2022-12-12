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

sample_7_MOI_500_data = Read10X(data.dir = "/raw_data/cellranger/count/7_MOI_500_GEX/outs/filtered_feature_bc_matrix")
sample_7_MOI_500 = CreateSeuratObject(counts = sample_7_MOI_500_data, project = "Sample_7_MOI_500", min.cells = 3, min.features = 200)
sample_7_MOI_500[["percent.mt"]] <- PercentageFeatureSet(sample_7_MOI_500, pattern = "^MT-")

sample_6_MOI_100_data = Read10X(data.dir = "/raw_data/cellranger/count/6_MOI_100_GEX/outs/filtered_feature_bc_matrix")
sample_6_MOI_100 = CreateSeuratObject(counts = sample_6_MOI_100_data, project = "Sample_6_MOI_100", min.cells = 3, min.features = 200)
sample_6_MOI_100[["percent.mt"]] <- PercentageFeatureSet(sample_6_MOI_100, pattern = "^MT-")

sample_5_HCT_116_data = Read10X(data.dir = "/raw_data/cellranger/count/5_HCT_116_GEX/outs/filtered_feature_bc_matrix")
sample_5_HCT_116 = CreateSeuratObject(counts = sample_5_HCT_116_data, project = "Sample_5_HCT_116", min.cells = 3, min.features = 200)
sample_5_HCT_116[["percent.mt"]] <- PercentageFeatureSet(sample_5_HCT_116, pattern = "^MT-")

sample.combine <- merge(sample_5_HCT_116, y = c(sample_6_MOI_100,sample_7_MOI_500), add.cell.ids = c("5_HCT_116_GEX","6_MOI_100_GEX",'7_MOI_500_GEX'), project = "SAMPLE.INTEGRATED")

VlnPlot(sample.combine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sample.combine <- NormalizeData(object = sample.combine, normalization.method = "LogNormalize", scale.factor = 10000)
GetAssay(sample.combine,assay = "RNA")
sample.combine <- FindVariableFeatures(object = sample.combine, selection.method = "vst", nfeatures = 5000)
top20 <- head(x = VariableFeatures(object = sample.combine), 20)
plot1 <- VariableFeaturePlot(object = sample.combine)
plot2 <- LabelPoints(plot = plot2, points = top20, repel = TRUE)
plot1+plot2
all.genes <- rownames(sample.combine)
sample.combine<- ScaleData(object = sample.combine,features = all.genes)
sample.combine <- RunPCA(object = sample.combine,pc.genes = VariableFeatures(sample.combine))
ElbowPlot(sample.combine)
seuratObj <- RunHarmony(sample.combine, group.by.vars="orig.ident",assay.use='RNA')
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:20, 
                     reduction = "harmony",seed.use=111)
sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:20)
sce <- FindClusters(sce, resolution = 0.5)
seuratObj=sce
sample.combine = seuratObj

# add pathogen UMI metadata
umi_table_csv = 'csv_novami_mix_dedup.csv'
umi_table<-read.csv(umi_table_csv,sep=',',header=TRUE,row.names = 1)
umi_table[is.na(umi_table)] <- 0
umi_table$Total <- rowSums(umi_table)
umi_table[umi_table==0] <- NA
sample.headneck<-AddMetaData(sample.headneck, umi_table)

