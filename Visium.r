# Load packages
library(ggplot2)
library(Seurat)
library(dplyr)
library(readr)
library(plyr)

setwd("/data_processing/R_processed_for_shiny")
path <- '/raw_data/count'
sample_pattern = 'V'# extract visium samples in the SpaceRanger folder
samples <- sort(list.files(path, pattern=sample_pattern, full.names = TRUE))
sample.names <- basename(samples)
processed_foler = '/data_processing/python'
samples_folder = '/raw_data/count'

output_path = '/data_processing/R_processed_for_shiny/rds_folder_clustering'

# generating rds files for each sample
for (each_sample in sample.names){

    print(each_sample)
    data_path = paste0(samples_folder,'/',each_sample,'/outs')
    metadata_file = paste0(processed_foler,'/',each_sample,'.visium.raw_matrix.genus.csv')
    output_file_filtered = paste0(output_path,'/',each_sample,'.filtered_matrix.rds')
    output_file_raw = paste0(output_path,'/',each_sample,'.raw_matrix.rds')

# processing: filtered matrix
    print('processing filtered matrix')
    list.files(samples_folder) # Should show filtered_feature_bc_matrix.h5
    tissue_sample<-Load10X_Spatial(data.dir = data_path,filename = "filtered_feature_bc_matrix.h5")
    plot1 <- VlnPlot(tissue_sample, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(tissue_sample, features = "nCount_Spatial") + theme(legend.position = "right")
    plot1+plot2
    tissue_sample <- subset(tissue_sample, subset = nCount_Spatial > 3 & nFeature_Spatial > 3)
    SpatialFeaturePlot(tissue_sample, features = "nCount_Spatial") + theme(legend.position = "right")
    tissue_sample <- SCTransform(tissue_sample, assay = "Spatial")
# delete original matrix to reduce size
    tissue_sample@assays$Spatial=NULL
    CellsMeta = tissue_sample@meta.data 
    head(CellsMeta)
    umi_table_csv = metadata_file
    umi_table<-read.csv(umi_table_csv,sep=',',header=TRUE,row.names = 1)
    umi_table$Total <- rowSums(umi_table)
    umi_table[umi_table==0] <- NA
# turn 0 to NA
    tissue_sample<-AddMetaData(tissue_sample, umi_table)

    p1 = SpatialFeaturePlot(tissue_sample,features = c('Total'),pt.size.factor = 1.52) +
                ggtitle(paste(each_sample,' ','Total Pathogen'," nUMI Filtered", sep = "")) + 
                theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
    print(p1)
# clustering

    tissue_sample <- RunPCA(tissue_sample, assay = "SCT", verbose = FALSE)
    tissue_sample <- FindNeighbors(tissue_sample, reduction = "pca", dims = 1:20)
    tissue_sample <- FindClusters(tissue_sample, verbose = FALSE)
    set.seed(123)
    tissue_sample <- RunUMAP(tissue_sample, reduction = "pca", dims = 1:20) 
    # find Spatially Variable Features
    tissue_sample <- FindSpatiallyVariableFeatures(tissue_sample, assay = "SCT", features = VariableFeatures(tissue_sample)[1:1000],selection.method = "markvariogram")

# save filtered rds
    saveRDS(tissue_sample, file = output_file_filtered)

# processing: raw matrix (mapping only)
    print('processing raw matrix')
    tissue_sample_raw<-Load10X_Spatial(data.dir = data_path,filename = "raw_feature_bc_matrix.h5")
    tissue_sample_raw<-AddMetaData(tissue_sample_raw, umi_table)

    p2 = SpatialFeaturePlot(tissue_sample_raw,features = c('Total'),pt.size.factor = 1.52) +
                ggtitle(paste(each_sample,' ','Total Pathogen'," nUMI RAW", sep = "")) + 
                theme(legend.position = "right",plot.title = element_text(hjust = 0.5))

    print(p2)
    saveRDS(tissue_sample_raw, file = output_file_raw)

}
