#load libraries 
library(Seurat) 
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)
#read and prepare seurat object
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
EC@assays[["integrated"]] <- NULL
Idents(EC) <- "PATH2"
setwd("F:/1- ATLAS/Figure 2/CellChat/3- pathological")
PATH = subset(EC, idents = "Pathology")
remove(EC)
DefaultAssay(PATH) <- "RNA"
Idents(PATH) <- "ECclusters"
levels(PATH) <- c("Large artery", "Artery", "Arteriole","Angiogenic capillary", "Capillary", "Venule", "Vein", "Large vein","Mitochondrial", "EndoMT","Proliferating EndoMT", "Proliferating cell","Stem-to-EC","Proliferating Stem-to-EC")
PATH$ECclusters <- PATH@active.ident


#-----------------------------------------------------------------------------#
#create Cellchat object from Seurat object
data.input <- GetAssayData(PATH, assay = "RNA", slot = "data") 
labels <- Idents(PATH)
meta <- data.frame(group = labels, row.names = names(labels)) 

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

remove(PATH)
remove(meta)
remove(labels)
# load the saved CellChatDB
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB  

#-----------------------------------------------------------------------------#
# CellChat analysis:
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#-----------------------------------------------------------------------------#
# Evalutuin of the interactions
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd','Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d','Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4','Proliferating Stem-to-EC' = '#89e56b')

#to plot circleplot.
groupSize <- as.numeric(table(cellchat@idents))

pdf("1- Number_of_interactions_circleplot.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = colours)
dev.off()
pdf("2- Interaction_weights_circleplot.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength", color.use = colours)
dev.off()
#-----------------------------------------------------------------------------#
# visualization of weights/strength of communications per cell type
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

#-----------------------------------------------------------------------------#
# MHC-II visualzation
#Visualization of cell-cell communication network -  MHC-II signaling pathway
pathways.show <- c("MHC-II")

pdf("5- MHC-II_circleplot.pdf", width = 8, height = 8)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = colours)
dev.off()

pdf("6- MHC-II_chordplot.pdf", width = 8, height = 8)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", color.use = colours)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, color.use = colours)
dev.off()

pdf("7- MHC-II_chordplot_Angiogenic_senders.pdf", width = 8, height = 8)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = FALSE, sources.use = 4, targets.use = c(1:14), color.use = colours)
dev.off()

pdf("8- MHC-II_chordplot_Large_veins_receivers.pdf", width = 8, height = 8)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", scale = TRUE, sources.use = c(1:14), targets.use = 8, color.use = colours)
dev.off()

pdf("9- MHC-II_heatmap.pdf", width = 8, height = 8)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", measure = "weight", slot.name = "netP", color.use = colours)
dev.off()
pdf("10- contribution.pdf", width = 8, height = 10)
  netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
pathways.show.all <- cellchat@netP$pathways
pdf("11- plotGeneExpression.pdf", width = 16, height = 12)
plotGeneExpression(cellchat, signaling = "MHC-II", color.use = colours)
dev.off()

#-----------------------------------------------------------------------------#
#Systems analysis of cell-cell communication network
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("12- Signaling_role_heatmap.pdf", width = 8, height = 4)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10, color.use = colours)
dev.off()

pdf("13- Signaling_role_scatter.pdf", width = 8, height = 8)
gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use = colours)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-II"), color.use = colours)
print(gg1 + gg2)
dev.off()


pdf("12- Outgoing_incoming_signaling_heatmaps.pdf", width = 12, height = 17)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 17, font.size = 5, color.use = colours)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 17, font.size = 5, color.use = colours)
print(ht1 + ht2)
dev.off()

#-----------------------------------------------------------------------------#

#Identify global communication patterns 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns, 
                                           heatmap.show = FALSE)
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming",
                                          k = nPatterns, heatmap.show = FALSE)


pdf("14 - River_plots.pdf", width = 8, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing", color.use = colours)
netAnalysis_river(cellchat, pattern = "incoming", color.use = colours)
dev.off()

pdf("15- Dot_plots.pdf", width = 12, height = 8)
netAnalysis_dot(cellchat, pattern = "outgoing", color.use = colours)
netAnalysis_dot(cellchat, pattern = "incoming", color.use = colours)
dev.off()

#-----------------------------------------------------------------------------#
# Functional and Structural 
# Functional
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")

cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE, nCores = 1)
pdf("16- Functional_embedding.pdf", width = 8, height = 8)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()
pdf("17- Functional_embedding_zoom.pdf", width = 8, height = 8)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()

# structural
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE, nCores = 1)
pdf("18- Structural_embedding.pdf", width = 8, height = 8)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf("19- Structural_embedding_zoom.pdf", width = 8, height = 8)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

#Save the CellChat object
saveRDS(cellchat, file = "Overall merge of pathological sorted endothelial cells_CellChat object.rds")
