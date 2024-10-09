library(CellChat)
library(Seurat)
#Prepare cellchat objects
TL <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/2- TL/Adult control brain (temporal lobe) sorted endothelial cells_CellChat object.rds")
PATH <- readRDS(file = "F:/1- ATLAS/Figure 2/CellChat/3- pathological/Overall merge of pathological sorted endothelial cells_CellChat object.rds")
setwd("F:/1- ATLAS/Figure 2/CellChat/4- control pathological")
object.list <- list(TL = TL, PATH = PATH)
remove(TL)
remove(PATH)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#------------------------------------------------------------------------#
#compare between them
pdf("1- compare interaction.pdf", width = 14, height = 8)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()
colours <- c('Large artery' = '#a1000b', 'Artery' = '#d0000d', 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd','Vein' = '#29e2e6', 'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d','Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4','Proliferating Stem-to-EC' = '#89e56b')

pdf("2- interactions number.pdf", width = 10, height = 8)

# Set up the plotting layout
netVisual_diffInteraction(cellchat, weight.scale = TRUE, color.use = colours)
dev.off()
pdf("3- interaction strengths")
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", color.use = colours)
dev.off()

pdf("4- netVisual_heatmap_combined.pdf", width = 10, height = 8)
gg1 <- netVisual_heatmap(cellchat, color.use = colours)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = colours)
gg1 + gg2
dev.off()

# For number of interactions
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

for (i in 1:length(object.list)) {
  pdf(file = paste0("5- Number_of_interactions_", names(object.list)[i], ".pdf"),width = 10,height = 10 )  # Open a PDF file
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, color.use = colours, 
                   label.edge= FALSE, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
  dev.off() 
}

# For strength of interactions
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))

for (i in 1:length(object.list)) {
  pdf(file = paste0("6- Strength_of_interactions_", names(object.list)[i], ".pdf"), width = 10,height = 10)  # Open a PDF file
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = TRUE, color.use = colours, 
                   label.edge= FALSE, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Strength of interactions - ", names(object.list)[i]))
  dev.off()  
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Control dot size in different datasets
pdf("7- outgoing_TL_PATH.pdf", height = 8, width = 14)
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax, 
                                               color.use = colours)
}
patchwork::wrap_plots(plots = gg)

# Close the PDF device after saving all plots
dev.off()
#------------------------------------------------------------------------#
# Define a vector of cell types
cell_types <- c("Large artery", "Artery", "Arteriole", "Capillary", "Venule", 
                "Vein", "Large vein", "Mitochondrial", "EndoMT", 
                "Proliferating EndoMT", "Proliferating cell", 
                "Stem-to-EC", "Proliferating Stem-to-EC")

# Loop through each cell type and create a scatter plot
for (cell_type in cell_types) {
  pdf(paste0("8-signalingChanges_scatter_", gsub(" ", "_", cell_type), ".pdf"), 
      width = 10, height = 8)
  
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = cell_type, 
                                       font.size = 15, label.size = 3, 
                                       font.size.title = 15)
  
  dev.off()
}

pdf("8-signalingChanges_scatter_Angiogenic_capillary.pdf", width = 10, height = 8)
netAnalysis_signalingChanges_scatter(cellchat, 
                                     idents.use = "Angiogenic capillary", 
                                     font.size = 15, 
                                     label.size = 3, font.size.title = 15)+
  theme(
    legend.position = c(0.65, 0.3),       # Position inside the plot (values from 0 to 1 for x and y axis)
    legend.justification = c(0, 1),      # Anchor the legend at a specific point (bottom-left inside)
    legend.background = element_rect(fill = "white", color = "white"), # Add background to the legend
    legend.box.background = element_rect(color = "white", size = 0.5)  # Border around the legend box
  )

dev.off()
#------------------------------------------------------------------------#
#Identify the conserved and context-specific signaling pathways
# functional
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE, nCores = 1)
pdf("9-embeddingpairwise_functional.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 1)
dev.off()

# Structural
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE, nCores = 1)
pdf("10- embeddingpairwise_structural.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 1)
dev.off()

#------------------------------------------------------------------------#
# Diffferential Pathways
pdf("11- rank_similarity.pdf")
rankSimilarity(cellchat, type = "functional")
dev.off()
col = c("TL"="turquoise3","PATH"="#f8766d")
pdf("12- ranknet.pdf", width = 10, height = 13)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use = col)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = col)
gg1 + gg2

dev.off()

#------------------------------------------------------------------------
library(ComplexHeatmap)
# outgoing signalling patterns
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 20, color.use = colours)
pdf("13- outgoing_signaling_patterns.pdf", width = 15, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# incoming signalling patterns
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 20, color.heatmap = "GnBu", color.use = colours)

pdf("14- incoming_signaling_patterns.pdf", width = 15, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

dev.off()

# overall signalling patterns
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20, color.heatmap = "OrRd", color.use = colours)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 20, color.heatmap = "OrRd", color.use = colours)
pdf("15- overall_signalling_patterns_heatmap.pdf", width = 16, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



#------------------------------------------------------------------------#
#compare communication probabilities in pathological (PATH) vs adult control brain (TL) endothelial cells
pdf("16- comm_prob_4.pdf", height = 40, width = 14)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14), 
                 comparison = c(1, 2), angle.x = 45)
dev.off()

pdf("17- comm_prob_1-14.pdf", height = 40, width = 14)
netVisual_bubble(cellchat, sources.use = c(1:14),
                 targets.use = 4,  comparison = c(1, 2), angle.x = 45)
dev.off()

# identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs 
pdf("18 - upregulated pathways.pdf", height = 40, width = 20)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14),
                 comparison = c(1, 2), max.dataset = 2,
                 title.name = "Increased signaling in PATH", 
                 angle.x = 45, remove.isolate = T) 
dev.off()
#> Comparing communications on a merged object
pdf("18 - downregulated pathways.pdf", height = 40, width = 20)

gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:14),
                        comparison = c(1, 2), max.dataset = 1,
                        title.name = "Decreased signaling in PATH",
                        angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
dev.off()

pdf("19- compare_comm_up_path.pdf", height = 40, width = 20)
netVisual_bubble(cellchat, sources.use = c(1:14),
                        targets.use =4,  comparison = c(1, 2),
                 max.dataset = 2, title.name = "Increased signaling in PATH",
                 angle.x = 45, remove.isolate = T)
dev.off()
# Comparing communications on a merged object
pdf("19- compare_comm_down_path.pdf", height = 40, width = 20)

netVisual_bubble(cellchat, sources.use = c(1:14), targets.use = 4, 
                 comparison = c(1, 2), max.dataset = 1,
                 title.name = "Decreased signaling in PATH",
                 angle.x = 45, remove.isolate = T)
dev.off()
#------------------------------------------------------------------------#
# Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "PATH"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with up-regulated ligands in PATH
net.up <- subsetCommunication(cellchat, net = net, datasets = "PATH",ligand.logFC = 0.1, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with up-regulated ligands and up-regulated receptors in TL, i.e.,down-regulated in PATH
net.down <- subsetCommunication(cellchat, net = net, datasets = "PATH",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                        sources.use = 4, targets.use = c(1:14),
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = 4, targets.use = c(1:14), 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pdf("20_upregulated.pdf", width = 15, height = 24)
gg1
dev.off()
pdf("21- downregulated.pdf", width = 15, height = 24)
gg2
dev.off()

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up
                        ,sources.use = c(1:14), targets.use = 4, 
                        comparison = c(1, 2), angle.x = 90, 
                        remove.isolate = TRUE,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = FALSE]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:14), targets.use = 4, comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


# Save the combined plot as a PDF
pdf("22- down-up in path.pdf", width = 8, height = 18)
gg1
dev.off()
pdf("23- upreguulated in path.pdf", width = 8, height = 10)
gg2
dev.off()
# Save first set of plots
pdf("24 - up_regulated_signaling_1.pdf", width = 15, height = 15)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, 
                     targets.use = c(1:14), slot.name = 'net',
                     color.use = colours, net = net.up, 
                     lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list)[2]))
netVisual_chord_gene(object.list[[2]], sources.use = c(1:14), targets.use = 4, slot.name = 'net',color.use = colours, net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

# Save plots for MHC-II pathway
pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
pdf("25- MHC_II_signaling.pdf", width = 20, height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# Repeat similar structure for VEGF, NOTCH, ANGPT, ANGPTL, APELIN pathways
pathways_list <- c("VEGF", "NOTCH", "ANGPT", "ANGPTL", "APELIN")
for (pathway in pathways_list) {
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathway)
  pdf(paste0(26, "_", pathway, "_signaling.pdf"), width = 20, height = 10)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathway, layout = "circle", color.use = colours, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
}

# heatmap for MHC-II pathway
pdf("27- MHC_II_heatmap.pdf")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], color.use = colours, signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# Save aggregated plots for MHC-II signaling
pdf("28- MHC_II_aggregated_signaling.pdf")
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, color.use = colours, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
# Save chord plots for differential analysis
sources <- c(3, 4, 5, 6, 7, 8, 10, 13)  # Adjust this vector as necessary
source_names <- c("Arteriole", "Angiogenic", "Capillary", "Venule", "Vein", "LargeVein", "ENDOMT", "Stem-to-EC")

for (i in seq_along(sources)) {
  pdf(paste0("29- MHC_II_chord_source_", source_names[i], ".pdf"), width = 14)
  
  netVisual_chord_gene(object.list[[2]], 
                       sources.use = sources[i],
                       color.use = colours,
                       signaling = pathways.show, 
                       targets.use = c(1:14), slot.name = 'net',
                       net = net.up, lab.cex = 0.8, 
                       small.gap = 3.5, title.name = paste0("MHC-II Up-regulated signaling in ", names(object.list)[2]))
  
  dev.off()
}

# EndoMT as target
pdf("30- MHC_II_chord_EndoMT_target.pdf", width = 14)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:14), color.use = colours,
                     signaling = pathways.show, targets.use = 10, slot.name = 'net',
                     net = net.up, lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("MHC-II Up-regulated signaling in ", names(object.list)[2], " (EndoMT)"))
dev.off()


# Save CellChat object
saveRDS(cellchat, file = "Pathological vs adult control brain sorted endothelial cells comparison_CellChat object.rds")


