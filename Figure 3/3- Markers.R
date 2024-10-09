library("heatmaply")
library(Seurat)
library(dplyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")

# set default assay to RNA
DefaultAssay(EC) <- "RNA"
EC@assays[["integrated"]] <- NULL
Idents(EC) <- "PATH3"

# [1] Capillary                Angiogenic capillary     Artery                  
# [4] Arteriole                Venule                   Vein                    
# [7] Mitochondrial            Stem-to-EC               Large vein              
# [10] Proliferating cell       Large artery             EndoMT                  
# [13] Proliferating Stem-to-EC Proliferating EndoMT 
# table(Capillary@meta.data$PATH3)

# MAL        TUM    CONTROL FETALbrain 
# Subset for malformations (MAL)
# Subset for Capillaries and each condition
# Function to find and save markers for specified conditions and vessel types
find_and_save_markers <- function(EC, condition, output_dir) {
  # Subset the data for the given condition
  subset_data <- subset(EC, subset = PATH3 == condition)
  Idents(subset_data) <- 'ECclusters'
  
  # Define vessel types
  vessel_types <- c("Capillary", "Large vein", "Large artery")
  
  # Loop through vessel types and find markers
  for (vessel in vessel_types) {
    markers <- FindMarkers(subset_data, ident.1 = vessel)
    # Create a file name for saving
    file_name <- paste("markers", tolower(gsub(" ", "_", vessel)), tolower(condition), ".csv", sep = "_")
    write.csv(markers, file = file.path(output_dir, file_name))
  }
}

# Set your output directory
output_directory <- "F:/1- ATLAS/Figure 3/M-O"

# Run the function for each condition
find_and_save_markers(EC, "MAL", output_directory)
find_and_save_markers(EC, "TUM", output_directory)
find_and_save_markers(EC, "CONTROL", output_directory)
find_and_save_markers(EC, "FETALbrain", output_directory)


#-----------------------------------------------------------------------------#

markers_capillary_mal$tissue_type <- rep("MAL", 
                                         length(rownames(markers_capillary_mal)))
markers_capillary_fetal$tissue_type <- rep("fetal", 
                                           length(rownames(markers_capillary_fetal)))

markers_Large_artery_mal$tissue_type <- rep("MAL", 
                                            length(rownames(markers_Large_artery_mal)))
markers_Large_artery_fetal$tissue_type <- rep("fetal", 
                                              length(rownames(markers_Large_artery_fetal)))

markers_Large_vein_mal$tissue_type <- rep("MAL", 
                                          length(rownames(markers_Large_vein_mal)))
markers_Large_vein_fetal$tissue_type <- rep("fetal", 
                                            length(rownames(markers_Large_vein_fetal)))


features_large_vein <- unique(c(markers_Large_vein_mal$X,markers_Large_vein_fetal$X ))
features_large_artery <- unique(c(markers_Large_artery_mal$X, markers_Large_artery_fetal$X))
features_capillary <- unique(c(markers_capillary_mal$X, markers_capillary_mal$X))


Idents(EC) <- "PATH3"
FETAL_MAL <- subset(EC, idents = c("FETALbrain","MAL"))
gene_expression <- GetAssayData(FETAL_MAL)
gene_expression_LV <- gene_expression[features_large_vein, ] 
gene_expression_LA <- gene_expression[features_large_artery,]
gene_expression_CAP <- gene_expression[features_capillary,]
occurrence_LV <- apply(gene_expression_LV, 1, function(x) {
  mean(x > 0) * 100  
})
occurrence_LA <- apply(gene_expression_LA, 1, function(x) {
  mean(x > 0) * 100  
})
occurrence_LA <- apply(gene_expression_LA, 1, function(x) {
  mean(x > 0) * 100  
})

