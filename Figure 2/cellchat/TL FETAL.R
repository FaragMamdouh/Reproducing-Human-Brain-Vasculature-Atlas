library(CellChat)
library(Seurat)
library(ggplot2)
#----------------------------------------------------------------------------#
# prepare the objects
TL <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/2- TL/Adult control brain (temporal lobe) sorted endothelial cells_CellChat object.rds")
FETAL <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/1- Fetal/Fetal_CNS_sorted_endothelial_cells_CellChat_object.rds")
setwd("F:/1- ATLAS/Figure 2/CellChat/5- Fetal Control")
object.list <- list(TL = TL, FETAL = FETAL)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#----------------------------------------------------------------------------#
# compare between objects
pdf("1- compareInteractions.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
print(gg1 + gg2)
dev.off()
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd','Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d','Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4','Proliferating Stem-to-EC' = '#89e56b')

pdf("2- diff_interaction_number.pdf", width = 12, height = 8)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = colours)
dev.off()

pdf("3- diff_interaction_strength.pdf", width = 12, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T,
                          measure = "weight", color.use = colours)
dev.off()

pdf("4- differential_heatmap.pdf", width = 14, height = 10)
gg1 <- netVisual_heatmap(cellchat,color.use = colours)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = colours)
gg1 + gg2
dev.off()
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("5- no of interaction featal TL.pdf", width = 16, height = 12)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, color.use = colours,label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("6- strength of interaction featal TL.pdf", width = 16, height = 12)
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T,color.use = colours,label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}
dev.off()



#----------------------------------------------------------------------------#
# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
pdf("7- source and outgoing fetal TL.pdf", width = 12, height = 8)
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = colours)
}
patchwork::wrap_plots(plots = gg)
dev.off()

#----------------------------------------------------------------------------#
# List of cell types to generate scatter plots for
cell_types <- c("Large artery", "Artery", "Arteriole", "Angiogenic capillary", 
                "Capillary", "Venule", "Vein", "Large vein", "Mitochondrial", 
                "EndoMT", "Proliferating EndoMT", "Proliferating cell", 
                "Stem-to-EC", "Proliferating Stem-to-EC")

start_figure_number <- 8

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  pdf(paste0( start_figure_number + i - 1, "_", cell_type, "_signalingChanges_scatter.pdf"))  
  p <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = cell_type, 
                                            font.size = 15, label.size = 3, font.size.title = 15)
  print(p)  
  dev.off()  
}
#----------------------------------------------------------------------------#
#Identify the conserved and context-specific signaling pathways
# Functional
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE, nCores = 1)
pdf("22- functional TL FETAL.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 1)
dev.off()
#----------------------------------------------------------------------------#
#Structural
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE, nCores = 1)
pdf("23- structural TL FETAL.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 1)
dev.off()

#----------------------------------------------------------------------------#
#Diff Pathways
pdf(paste0(24, "_rankSimilarity.pdf"))
rankSimilarity(cellchat, type = "functional")
dev.off()

col <- c("TL" = "turquoise3", "FETAL" = "#A0E4A8")
pdf(paste0(25, "_stacked_rankNet.pdf"), height = 12)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, color.use = col)
print(gg1)
dev.off()

pdf(paste0(26, "_unstacked_rankNet.pdf"), height = 12)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, color.use = col)
print(gg2)
dev.off()

pdf(paste0(27, "_combined_rankNet.pdf"), height = 12)
gg1 + gg2
dev.off()
#----------------------------------------------------------------------------#
# Outgoing and incoming
# Plot heatmap of signaling patterns - outgoing
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)

pdf(paste0(28, "_outgoing_signaling_heatmap.pdf"), height = 10, width = 10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.use = colours)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(paste0(29, "_incoming_signaling_heatmap.pdf"), height = 10, width = 10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(paste0(30, "_all_signaling_heatmap.pdf"), height = 10, width = 10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "OrRd", color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i + 1], width = 10, height = 20, color.heatmap = "OrRd", color.use = colours)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# Compare communication probabilities in fetal vs adult control brain
pdf(paste0(31, "_bubble_comparison_fetal_vs_adult.pdf"), height = 40, width = 8)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14), 
                 comparison = c(1, 2), angle.x = 45) 
dev.off()

pdf(paste0(32, "_bubble_comparison_adult_vs_fetal.pdf"),height = 40, width = 8 )
netVisual_bubble(cellchat, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), angle.x = 45)
dev.off()

# Up-regulated and down-regulated signaling ligand-receptor pairs
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PATH", angle.x = 45, remove.isolate = TRUE)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PATH", angle.x = 45, remove.isolate = TRUE)

pdf(paste0(33, "_up_down_regulated_signaling_1.pdf"), height = 40, width = 16)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PATH", angle.x = 45, remove.isolate = TRUE)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PATH", angle.x = 45, remove.isolate = TRUE)

pdf(paste0(34, "_up_down_regulated_signaling_2.pdf"), height = 30, width = 12)
gg1 + gg2
dev.off()
#----------------------------------------------------------------------------#
# Identify dysfunctional signaling
pos.dataset = "FETAL"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)

# Extract up-regulated ligand-receptor pairs
net.up <- subsetCommunication(cellchat, net = net, datasets = "FETAL", ligand.logFC = 0.1, receptor.logFC = 0.1)
# Extract down-regulated ligand-receptor pairs
net.down <- subsetCommunication(cellchat, net = net, datasets = "FETAL", ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = FALSE]
pdf(paste0(35, "_up_regulated_signaling_1.pdf"), height = 18)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(1:14), comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
print(gg1)
dev.off()

pairLR.use.down = net.down[, "interaction_name", drop = FALSE]
pdf(paste0(36, "_down_regulated_signaling_1.pdf"))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(1:14), comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
print(gg2)
dev.off()


pairLR.use.up = net.up[, "interaction_name", drop = FALSE]
pdf(paste0(38, "_up_regulated_signaling_2.pdf"), height = 18)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
print(gg1)
dev.off()

pairLR.use.down = net.down[, "interaction_name", drop = FALSE]
pdf(paste0(39, "_down_regulated_signaling_2.pdf"))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
print(gg2)
dev.off()

pdf(paste0(40, "_up_down_regulated_signaling_2_combined.pdf"), width = 22, height = 30)
gg1 + gg2
dev.off()


saveRDS(cellchat, file = "Fetal vs adult control brain sorted endothelial cells comparison_CellChat object.rds")

