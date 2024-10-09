library(patchwork)
#----------------------------------------------------------------------------#
# Prepare CellObjects
FETAL <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/1- Fetal/Fetal_CNS_sorted_endothelial_cells_CellChat_object.rds")
TL <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/2- TL/Adult control brain (temporal lobe) sorted endothelial cells_CellChat object.rds")
PATH <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/3- pathological/Overall merge of pathological sorted endothelial cells_CellChat object.rds")
setwd("F:/1- ATLAS/Figure 2/CellChat/6- All comparsion")
object.list <- list(FETAL = FETAL, TL = TL, PATH = PATH)
remove(FETAL)
remove(TL)
remove(PATH)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#----------------------------------------------------------------------------#
# Objects comparsion
# Set up colors for plots
col = c("FETAL"="#A0E4A8", "TL"="turquoise3", "PATH"="#f8766d")
pdf("1_cellchat_analysis_figures.pdf", width = 10, height = 8)

# Compare interactions and plot
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3), color.use = c("#A0E4A8", "turquoise3", "#f8766d"))  
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2, 3), measure = "weight", color.use = c("#A0E4A8", "turquoise3", "#f8766d")) 
print(gg1 + gg2)  
dev.off()

# Set up additional colors for endothelial cell types
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 
             'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd', 
             'Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 
             'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d', 
             'Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4', 
             'Proliferating Stem-to-EC' = '#89e56b')

# Get max weight for count
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(1, 3), xpd = TRUE)

# Loop through object.list to create plots for number of interactions
for (i in 1:length(object.list)) {
  pdf(file = paste0("2_Number_of_interactions_", names(object.list)[i], ".pdf"), width = 12, height = 10)  
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = TRUE, 
                   color.use = colours, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
  dev.off() 
}

# Get max weight for weight and set up the plotting area
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
par(mfrow = c(1, 3), xpd = TRUE)

# Loop through object.list to create plots for strength of interactions
for (i in 1:length(object.list)) {
  pdf(file = paste0("3_Strength_of_interactions_", names(object.list)[i], ".pdf"), width = 12, height = 10)  # Open PDF device
  netVisual_circle(object.list[[i]]@net$weight, 
                   weight.scale = TRUE, 
                   color.use = colours, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Strength of interactions - ", names(object.list)[i]))
  dev.off() 
}
#----------------------------------------------------------------------------#
# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) { rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count) })
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()

# Create scatter plots for each object in object.list
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax, 
                                               color.use = colours)
}

pdf(file = "4_Combined_Scatter_Plots.pdf", width = 12, height = 8)  
wrap_plots(plots = gg)
dev.off()  
#----------------------------------------------------------------------------#
#Identify the conserved and context-specific signaling pathways
# Functional
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE, nCores = 1)
pdf("5- functional.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 2)
dev.off()
# Structural
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE, nCores = 1)
pdf("6- structural.pdf", width = 10)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 2)
dev.off()
#----------------------------------------------------------------------------#
library(ComplexHeatmap)
# Plot heatmap of signaling patterns
i = 1
# Get union of pathways from the first three objects
pathway.unionx <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
pathway.union <- union(pathway.unionx, object.list[[i + 2]]@netP$pathways)

# Plot outgoing signaling heatmaps
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing",
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], width = 10, 
                                        height = 20, color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
                                        pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i + 1], 
                                        width = 10, height = 20,
                                        color.use = colours)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i + 2]], pattern = "outgoing", 
                                        signaling = pathway.union,
                                        title = names(object.list)[i + 2], width = 10,
                                        height = 20, color.use = colours)
pdf(file = "7_Outgoing_Signaling_Heatmap.pdf", height = 20, width = 30)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.off()

# Plot incoming signaling heatmaps
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i + 2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 2], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
pdf(file = "8_Incoming_Signaling_Heatmap.pdf", width = 20, height = 30)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.off()



#----------------------------------------------------------------------------#
# Plot all signaling heatmaps
ht1_FETAL = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                              pattern = "all",
                                              signaling = pathway.union,
                                              title = names(object.list)[i],
                                              color.heatmap = "OrRd",
                                              color.use = colours, width = 15, 
                                              height = 30)

png("9_All_Signaling_Heatmap_FETAL.png", width = 800,height = 1600)
ht1_FETAL
dev.off()
ht2_TL = netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
                                           pattern = "all", signaling = pathway.union,
                                           title = names(object.list)[i + 1], 
                                           width = 15, height = 30, 
                                           color.heatmap = "OrRd", color.use = colours)+
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))
png("9_All_Signaling_Heatmap_TL.png", width = 800, height = 1600)
ht2_TL
dev.off()
ht3_PATH = netAnalysis_signalingRole_heatmap(object.list[[i + 2]], pattern = "all",
                                             signaling = pathway.union, 
                                             title = names(object.list)[i + 2],
                                             color.heatmap = "OrRd",
                                             color.use = colours, width = 15, 
                                             height = 30)+
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))
png("9_All_Signaling_Heatmap_PATH.png", width = 800, height = 1600)
ht3_PATH
dev.off()

#----------------------------------------------------------------------------#
# Visually compare cell-cell communication using Hierarchy plot for MHC-II pathway
pathways.show <- c("MHC-II")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

# Loop through object.list to create and save MHC-II pathway plots
for (i in 1:length(object.list)) {
  pdf(file = paste0("10_MHC-II_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off()  
}

# Repeat for VEGF pathway
pathways.show <- c("VEGF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("11_VEGF_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off() 
}

# Repeat for NOTCH pathway
pathways.show <- c("NOTCH")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("12_NOTCH_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off() 
}

# Repeat for ANGPT pathway
pathways.show <- c("ANGPT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("13_ANGPT_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off()  
}

# Repeat for ANGPTL pathway
pathways.show <- c("ANGPTL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("14_ANGPTL_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off()  
}

# Repeat for APELIN pathway
pathways.show <- c("APELIN")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Control the edge weights across different datasets
par(mfrow = c(1, 3), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("15_APELIN_Pathway_Aggregate_", names(object.list)[i], ".pdf"), width = 9, height = 9)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off() 
}

# Heatmap for MHC-II pathway
pathways.show <- c("MHC-II")
par(mfrow = c(1, 3), xpd = TRUE)
ht <- list()

for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], color.use = colours, signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(object.list)[i]))
}

# Save heatmap for MHC-II pathway
pdf(file = "16_MHC-II_Pathway_Heatmap.pdf", width = 12)
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))
dev.off()

# Aggregate for MHC-II pathway in chord diagram
par(mfrow = c(1, 1), xpd = TRUE)

for (i in 1:length(object.list)) {
  pdf(file = paste0("17_MHC-II_Pathway_Chord_", names(object.list)[i], ".pdf"), width = 8, height = 6)  # Open PDF device
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, color.use = colours, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
  dev.off()  
}

saveRDS(cellchat, file = "Fetal vs adult control vs Pathological brain sorted endothelial cells comparison_CellChat object.rds")

