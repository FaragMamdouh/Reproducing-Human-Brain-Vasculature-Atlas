# Load libraries 
library(Seurat) 
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
# Set working directory and load Seurat object
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
setwd("F:/1- ATLAS/Figure 2/CellChat/1- Fetal")
Idents(EC) <- "PATH"
FETAL <- subset(EC, idents = "FETALCNS")
remove(EC)
DefaultAssay(FETAL) <- "RNA"
Idents(FETAL) <- "ECclusters"
levels(FETAL) <- c("Large artery", "Artery", "Arteriole","Angiogenic capillary", "Capillary", "Venule", "Vein", "Large vein","Mitochondrial", "EndoMT","Proliferating EndoMT", "Proliferating cell","Stem-to-EC","Proliferating Stem-to-EC")
FETAL$ECclusters <- FETAL@active.ident
#------------------------------------------------------------------------------#

#create Cellchat object from Seurat object
data.input <- GetAssayData(FETAL, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(FETAL)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

remove(FETAL)
remove(data.input)
remove(meta)
#load the saved CellChat database
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB  
remove(CellChatDB)
#------------------------------------------------------------------------------#

# Preprocessing and computing communication network
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#------------------------------------------------------------------------------#
# Visualize the results
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd','Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d','Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4','Proliferating Stem-to-EC' = '#89e56b')

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("1- Interaction_numbder_circleplot.pdf",width = 10, height = 10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = colours)
dev.off()

pdf("2- Interaction_strength_circleplot.pdf", width = 10, height = 10)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = colours)
dev.off()

#------------------------------------------------------------------------------#
# per cell type interactions

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0("3- Interaction_strength_", rownames(mat)[i], ".pdf"), width = 10, height = 10)
  netVisual_circle(mat2, vertex.weight = groupSize, color.use = colours, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

mat <- cellchat@net$count
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0("4- Interaction_count_", rownames(mat)[i], ".pdf"), width = 10, height = 10)
  netVisual_circle(mat2, vertex.weight = groupSize, color.use = colours, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}


#------------------------------------------------------------------------------#
# MHC-II signaling pathway visualizations
pathways.show <- c("MHC-II")

pdf("5- MHC-II_circleplot.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = colours)
dev.off()

pdf("6- MHC-II_chordplot.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = colours)
dev.off()

pdf("7- MHC-II_chordplot_scaled.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, color.use = colours)
dev.off()

pdf("8- MHC-II_chordplot_sender_AngiogenicCapillary.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = FALSE, sources.use = 4, targets.use = c(1:14), color.use = colours)
dev.off()

pdf("9- MHC-II_chordplot_receiver_LargeVeins.pdf")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, sources.use = c(1:14), targets.use = 8, color.use = colours)
dev.off()

pdf("10 -MHC-II_heatmap.pdf")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", measure = "weight", slot.name = "netP", color.use = colours)
dev.off()

# Contribution of ligand-receptor pairs
pdf("11- MHC-II_contribution.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

# Violin plot of signaling gene expression
pdf("12- MHC-II_violinplot.pdf")
plotGeneExpression(cellchat, signaling = "MHC-II", color.use = colours)
dev.off()

#------------------------------------------------------------------------------#
# Outgoing and incoming
cellchat <- netAnalysis_computeCentrality(cellchat)
pdf("13- Outgoing_signalingRole_heatmap.pdf", width = 8, height = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 17, font.size = 5, color.use = colours)
print(ht1)
dev.off()

pdf("13- Incoming_signalingRole_heatmap.pdf", width = 8, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 17, font.size = 5, color.use = colours)
print(ht2)
dev.off()
library(ggalluvial)

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
nPatterns = 3
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 8, height = 14,font.size = 5)
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 6, height = 14,font.size = 5)

pdf("14- Outgoing_riverplot.pdf")
netAnalysis_river(cellchat, pattern = "outgoing", color.use = colours)
dev.off()

pdf("14- Incoming_riverplot.pdf")
netAnalysis_river(cellchat, pattern = "incoming", color.use = colours)
dev.off()

pdf("15- Incoming_dotplot.pdf", width = 12)
netAnalysis_dot(cellchat, pattern = "incoming", color.use = colours)
dev.off()
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")


# Saving the CellChat object
saveRDS(cellchat, file = "Fetal_CNS_sorted_endothelial_cells_CellChat_object.rds")
