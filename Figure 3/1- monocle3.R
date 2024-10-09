# Load libraries
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggbeeswarm)
library(viridis)

# A function for Monocle3 pseudotime analysis and plotting
run_monocle_analysis <- function(seurat_object_path, output_dir, output_prefix ) {
  # Load the Seurat object
  cds <- readRDS(file = seurat_object_path)
  
  # Set working directory
  setwd(output_dir)
  
  # Plot the pseudotime UMAP
  output_file <- paste0(output_prefix, "_pseudotime.png")  # Specify output file as PNG
  
  # Create the plot using plot_cells
  pseudotime_plot <- plot_cells(cds, 
                                color_cells_by = "pseudotime", 
                                label_cell_groups = FALSE, 
                                label_leaves = FALSE, 
                                label_branch_points = FALSE, 
                                graph_label_size = 6,
                                group_label_size = 6,
                                cell_size = 1.5) + 
    scale_color_viridis(discrete = FALSE, option = "E")
  
  # Save the plot using ggsave
  ggsave(output_file, plot = pseudotime_plot, width = 10, height = 8, dpi = 300)
  
  # Plot cell order according to pseudotime
  pdata_cds <- pData(cds)
  pdata_cds$pseudotime_monocle3 <- monocle3::pseudotime(cds)
  
  # Define the colors for cell types and level order
  my_color <- c('EndoMT' = '#70fb52', 'Stem-to-EC' = '#4cd7a4', 
                'Mitochondrial' = '#f1b77c', 'Angiogenic capillary' = '#c300fe', 
                'Capillary' = '#fb00c9', 'Large artery' = '#a1000b', 'Artery' = '#d0000d',
                'Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Venule' = '#6dc4fd', 
                'Arteriole' = '#d13b18', 'Proliferating cell' = '#FBFF00', 
                'Proliferating Stem-to-EC' = '#89e56b', 'Proliferating EndoMT' = '#affd2d')
  
  level_order <- c('Large artery', 'Artery', 'Arteriole', 'Angiogenic capillary', 
                   'Capillary', 'Venule', 'Vein', 'Large vein', 'Mitochondrial', 
                   'EndoMT', 'Proliferating EndoMT', 'Proliferating cell', 
                   'Stem-to-EC', 'Proliferating Stem-to-EC')
  
  output_file <- paste0(output_prefix, "_cell_order.png")
  plot <- ggplot(as.data.frame(pdata_cds), 
                 aes(x = pseudotime_monocle3, 
                     y = factor(ECclusters, level = level_order), 
                     colour = ECclusters)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_manual(values = my_color) + 
    theme_classic() +
    xlab("pseudotime") + 
    ylab("CellType") 
  
  ggsave(output_file, plot = plot, width = 10, height = 8, dpi = 300)
    if(output_prefix == "TL"){
      ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
      write.csv(ciliated_cds_pr_test_res, "ciliated.csv")
      pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
      write.csv(pr_deg_ids, "pr_deg_ids.csv")
    }
}
# Visualize fetal condition
run_monocle_analysis(
  seurat_object_path = "C:/Users/Farag/Downloads/Fetal CNS sorted endothelial cells_cds monocle object.rds", 
  output_dir = "F:/1- ATLAS/Figure 3/1 monocle3/FETAL", 
  output_prefix = "fetal"
)
# Visualize TL condition
run_monocle_analysis(
  seurat_object_path = "C:/Users/Farag/Downloads/Adult control brain (temporal lobe) sorted endothelial cells_cds monocle object.rds", 
  output_dir = "F:/1- ATLAS/Figure 3/1 monocle3/TL", 
  output_prefix = "TL"
)
# Visualize Path condition
run_monocle_analysis(
  seurat_object_path = "C:/Users/Farag/Downloads/Overall merge of pathological sorted endothelial cells_cds monocle object.rds", 
  output_dir = "F:/1- ATLAS/Figure 3/1 monocle3/pathology", 
  output_prefix = "pathology"
)
