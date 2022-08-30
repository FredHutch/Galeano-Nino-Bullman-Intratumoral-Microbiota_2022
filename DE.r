#!/usr/bin/env Rscript
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")

m_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)


DE_GSEA <- function(seurat_object,
                                ident_1,
                                ident_2,
                                group_by,
                                seurat_object.markers_filename,
                                seurat_object.markers_filtered_filename,
                                seurat_object.markers_gsea_filename){
    seurat_object.markers <- FindMarkers(seurat_object, 
                                            ident.1 = ident_1,
                                            ident.2 = ident_2,
                                            group.by = group_by, 
                                            logfc.threshold = -Inf, 
                                            min.pct = 0.1)

    write.csv(seurat_object.markers,seurat_object.markers_filename, row.names = TRUE)

    #seurat_object.markers = filter(seurat_object.markers, p_val_adj <= 0.05)
    seurat_object.markers= seurat_object.markers[order(-seurat_object.markers$avg_log2FC),]
    seurat_object.markers_filename = seurat_object.markers_filtered_filename
    write.csv(seurat_object.markers,seurat_object.markers_filename, row.names = TRUE)

    markers_seurat_object <- seurat_object.markers[,c("avg_log2FC")]
    names(markers_seurat_object) = as.character(rownames(seurat_object.markers))
    markers_seurat_object = sort(markers_seurat_object, decreasing = TRUE)
    length(markers_seurat_object)

    markers_seurat_object.em2 <- GSEA(markers_seurat_object, 
                                        TERM2GENE = m_H,
                                        eps=0.0,
                                        by = "fgsea")

    write.csv(markers_seurat_object.em2,seurat_object.markers_gsea_filename, row.names = FALSE)
}

