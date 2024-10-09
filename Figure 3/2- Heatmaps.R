library("heatmaply")
library(Seurat)
library(dplyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
library("pheatmap")
# Set working directory and load data
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")

# Set default assay to RNA
DefaultAssay(EC) <- "RNA"
EC@assays[["integrated"]] <- NULL
Idents(EC) <- "PATH2"



# Define the function to generate heatmaps
generate_heatmap <- function(seurat_object, condition, features,order, output_dir) {
  # Subset the data based on condition
  subset_data <- subset(seurat_object, idents = condition)
  
  # Get the average expression per cluster and convert to data frame
  AV_markers <- AverageExpression(subset_data, assays = "RNA", features = features, 
                                  return.seurat = FALSE, group.by = "ECclusters",
                                  add.ident = NULL, slot = "data", verbose = TRUE)
  
  markers_dataframe <- as.data.frame(AV_markers)
  
  # Arrange the dataframe according to the desired order
  AV_markers_dataframe2 <- markers_dataframe[, order]
  
  colnames(AV_markers_dataframe2) <- gsub("^RNA\\.", "", colnames(AV_markers_dataframe2))
  
  # Generate heatmap and save as PDF
  output_file <- paste0(output_dir, "/AV_markers_Heatmap_", condition, ".pdf")
  pdf(output_file)
  pheatmap(AV_markers_dataframe2, 
           scale = "row",  
           cluster_rows = FALSE,  
           cluster_cols = FALSE,  
           show_colnames = TRUE,
           show_rownames = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
}

# Define features for heatmap
features <- rev(c(
  # Large Artery Markers
  "LTBP4", "ELN", "MGP", "BMX", "DKK2", "FBLN5",
  # Artery Markers
  "GLUL", "ADAMTS1", "ALPL", "GJA5", "VEGFC", "SEMA3G",
  # Arteriole Markers
  "AIF1L", "CD320", "NET1", "LGALS3", "VSIR",
  # Angiogenic Capillary Markers
  "PLVAP", "CA2", "ADM", "LXN", "PGF", "PXDN", "CXCR4", "ANGPT2", "APLN", "ESM1",
  # Capillary Markers
  "BSG", "SLC16A1", "SLCO1A2", "SLC38A5", "MFSD2A", "CA4",
  # Venule Markers
  "PRCP", "PRSS23", "RAMP3", "DNASE1L3", "POSTN", "PTGDS", "JAM2",
  # Vein Markers
  "NR2F2", "ACKR1", "CCL2", "IL1R1",
  # Large Vein Markers
  "CCL2", "VCAM1", "ICAM1", "SELP", "SELE"
))

# Output directory
output_dir <- "F:/1- ATLAS/Figure 3/D Heatmap"

condition <- "Control"
order <- c("RNA.Large.artery","RNA.Artery", 
           "RNA.Arteriole","RNA.Angiogenic.capillary", 
           "RNA.Capillary","RNA.Venule","RNA.Vein",
           "RNA.Large.vein","RNA.Mitochondrial", 
           "RNA.EndoMT",
           "RNA.Proliferating.cell",
           "RNA.Stem.to.EC")
#generate a heatmap for fetal condition
generate_heatmap(EC, condition, features, order, output_dir)
remove(Control)

order <- c("RNA.Large.artery","RNA.Artery", 
           "RNA.Arteriole","RNA.Angiogenic.capillary", 
           "RNA.Capillary","RNA.Venule","RNA.Vein",
           "RNA.Large.vein","RNA.Mitochondrial", 
           "RNA.EndoMT",
           "RNA.Proliferating.EndoMT",
           "RNA.Proliferating.cell",
           "RNA.Stem.to.EC",
           "RNA.Proliferating.Stem.to.EC")
#generate a heatmap for control and pathology condition
conditions  <- c("Fetal", "Pathology")
for (condition in conditions) {
  generate_heatmap(EC, condition, features,order, output_dir)
}
remove(Fetal)
remove(Pathology)
