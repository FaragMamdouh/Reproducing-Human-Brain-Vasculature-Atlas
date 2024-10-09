#load libraries 
library(Seurat) 
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
#read and prepra seurat object
options(stringsAsFactors = FALSE)
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
Idents(EC) <- "PATH"
TL <- subset(EC, idents = "TL")
setwd("F:/1- ATLAS/Figure 2/CellChat/2- TL")
remove(EC)
DefaultAssay(TL) <- "RNA"
Idents(TL) <- "ECclusters"
levels(TL) <- c("Large artery", "Artery", "Arteriole","Angiogenic capillary", "Capillary", "Venule", "Vein", "Large vein","Mitochondrial", "EndoMT", "Proliferating cell","Stem-to-EC")
TL$ECclusters <- TL@active.ident
#------------------------------------------------------------------------------#
#create Cellchat object from Seurat object
data.input <- GetAssayData(TL, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(TL)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat2 <- createCellChat(object = data.input, meta = meta, group.by = "group")
remove(TL)
remove(meta)
remove(data.input)

# load the saved CellChat database
CellChatDB <- CellChatDB.human 
cellchat2@DB <- CellChatDB  
remove(CellChatDB)
#------------------------------------------------------------------------------#

# CellChat analysis:
ellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2, type = "truncatedMean", trim = 0.05, population.size = TRUE)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)

#------------------------------------------------------------------------------#
#cellchatt visualzation
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd','Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d','Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4','Proliferating Stem-to-EC' = '#89e56b')

group = c("Large artery","Artery","Arteriole","Angiogenic capillary","Capillary","Venule","Vein","Large vein","Mitochondrial","EndoMT","Proliferating EndoMT","Proliferating cell",  "Stem-to-EC","Proliferating Stem-to-EC")
cellchat2 <- liftCellChat(cellchat2, group)

cellchat = cellchat2
# Plot circleplot
groupSize <- as.numeric(table(cellchat@idents))

pdf("1- circleplot_number_of_interactions.pdf", width = 10, height = 10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = colours)
dev.off()

pdf("2- circleplot_interaction_weights.pdf", width = 10, height = 10)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = colours)
dev.off()


#------------------------------------------------------------------------------#
#per cell type visualization
# To visualize weights/strength of communications per cell type
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0("3- Interaction_strength_", rownames(mat)[i], ".pdf"), width = 10, height = 10)
  netVisual_circle(mat2, vertex.weight = groupSize, color.use = colours, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
  
  }

# To visualize number of communications per cell type
mat <- cellchat@net$count
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0("4- Interaction_count_", rownames(mat)[i], ".pdf"), width = 10, height = 10)
  netVisual_circle(mat2, vertex.weight = groupSize, color.use = colours, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}


#------------------------------------------------------------------------------#
# MHC-II signaling pathway
pathways.show <- c("MHC-II")

pdf("5- MHC-II_circleplot.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = colours)
dev.off()

pdf("6- MHC-II_chordplot.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = colours)
dev.off()

pdf("7- MHC-II_chordplot_scaled.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, color.use = colours)
dev.off()

pdf("8- MHC-II_chordeplot_angiogenic_capillaries.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = FALSE, sources.use = 4, targets.use = c(1:14), color.use = colours)
dev.off()

pdf("9- MHC-II_chordeplot_large_veins.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, sources.use = c(1:14), targets.use = 8, color.use = colours)
dev.off()

pdf("10 -MHC-II_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", measure = "weight", slot.name = "netP", color.use = colours)
dev.off()

# Contribution of each ligand-receptor pair to the overall signaling pathway
pdf("11- netAnalysis_contribution.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pathways.show.all <- cellchat@netP$pathways

# Violin plot for signaling gene expression distribution
pdf("12 - MHC-II_violinplot.pdf")
plotGeneExpression(cellchat, signaling = "MHC-II", color.use = colours)
dev.off()
#------------------------------------------------------------------------------#
# Systems analysis of cell-cell communication network
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Centrality scores heatmap
pdf("13- signalingRole_network_heatmap.pdf")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10, color.use = colours)
dev.off()

# Scatter plot of dominant senders and receivers
pdf("14- signalingRole_scatter.pdf")
gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use = colours)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-II"), color.use = colours)
plot(gg1 + gg2)
dev.off()

# Heatmap of outgoing and incoming signaling roles
pdf("15- signalingRole_outgoing_heatmap.pdf", width = 8, height = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 17, font.size = 5, color.use = colours)
dev.off()

pdf("16- signalingRole_incoming_heatmap.pdf", height = 12, width = 8)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 17, font.size = 5, color.use = colours)
dev.off()

# Global communication patterns
colours1 <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd', 'Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4')

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing",
                                          k = nPatterns,
                                          heatmap.show = FALSE)

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = nPatterns, heatmap.show = FALSE)

# River plot
pdf("17- riverplot_outgoing.pdf")
netAnalysis_river(cellchat, pattern = "outgoing", color.use = colours)
dev.off()

pdf("18- riverplot_incoming.pdf")
netAnalysis_river(cellchat, pattern = "incoming", color.use = colours)
dev.off()

# Dot plot
pdf("19- dotplot_outgoing.pdf", width = 12)
netAnalysis_dot(cellchat, pattern = "outgoing", color.use = colours)
dev.off()

pdf("20- dotplot_incoming.pdf", width = 12)
netAnalysis_dot(cellchat, pattern = "incoming", color.use = colours)
dev.off()


#------------------------------------------------------------------------------#
# Functional and structural analysis
# Functional
pdf("21- functional_net_embedding.pdf")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional", 
                          do.parallel = FALSE, nCores = 1)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()
pdf("22- functional_zoomin.pdf")
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()


# Structural
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE, nCores = 1)
pdf("23- structural_net_embedding.pdf")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
pdf("24- structural_zoomin.pdf")
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()


#Save the CellChat object
saveRDS(cellchat, file = "Adult control brain (temporal lobe) sorted endothelial cells_CellChat object.rds")

