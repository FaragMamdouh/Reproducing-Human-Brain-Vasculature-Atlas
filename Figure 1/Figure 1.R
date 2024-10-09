
#Top entity (endothelial) markers heatmap - Figure 1b
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


#-------------------------------Figure 1b -------------------------------------#
# Input the the top entity markers features to plot
features = c("TTN","PRXL2A","TFRC" ,"NPIPB5","SLC7A1",
             "ADIRF","IFI27","LY6E","LGALS3", "CD320",
             "MGP","CCL2","SOD2", "TM4SF1","ACKR1",
             "CCN2","PLCG2","IFNGR1","EPAS1","TACC1", 
             "ANGPT2","CA2","SPINK8", "APLN","CTGF",
             "CYTOR","INHBB","INSR","SERPINA5", "ADAMTS5",
             "RBP7","GSN","TXNIP","PLPP1","AQP1")

markers <- AverageExpression(EC, assays = "RNA", features = features, return.seurat = FALSE, group.by = "PATH", add.ident = NULL, slot = "data", verbose = TRUE)
markers_dataframe <- as.data.frame(markers)

markers_dataframe2 <- as.data.frame(markers)[,c("RNA.FETALCNS",
                                                "RNA.TL","RNA.AVM","RNA.LGG",
                                                "RNA.GBM","RNA.MET","RNA.MEN")]
colnames(markers_dataframe2) <- c("Fetal Brain", "TL", "AVM", "LGG","GBM", "MET","MEN")
heatmaply(markers_dataframe2, scale = "row", colors = c("#0571b0","white","#ca0020"),
          grid_gap = 1,dendrogram = ("none"),
          plot_method = "plotly")


#-------------------------------Figures 1C-F-----------------------------------#
features1 = c("ATP5F1E", "PRXL2A", "TTN", "COL4A1", "SEPTIN7", "KDR", "PALM2-AKAP2",
              "MARCKS", "SEPTIN2", "CD93", "COL4A2", "NREP", "MYO1B", "SOX4", "RELL1",
              "NPIPB5", "APCDD1", "GTF2I", "SLC38A2", "MACF1", "IGFBP2", "CEMIP2",
              "LAMA4", "SMG1", "APLNR", "RNASE1", "NET1", "CAVIN2", "RPS17", "ATP5L",
              "ATP5I", "MT1M", "HERPUD1", "LY6E", "CTGF", "TSC22D1", "HLA-A", "SRGN",
              "LGALS3", "BST2", "B2M", "MT1E", "HLA-C", "MALAT1-ENSG00000251562",
              "ATP5E", "SPARCL1", "HLA-B", "MT2A", "ADIRF", "IFI27")

features2 = c("SPRY1", "MGP", "HSPG2", "TIMP1", "COL4A1", "PLVAP", "COL4A2", "CD93",
              "CALCRL", "SPARC", "ANGPT2", "STC1", "INSR", "VWA1", "PDLIM1", "MARCKSL1",
              "IGFBP4", "IGFBP3", "HTRA1", "EMP1", "VWF", "MCAM", "DUSP6", "IGFBP7",
              "ARHGDIB", "IFI6", "ABCG2", "TSC22D1", "COBLL1", "SRARP", "SLC38A5",
              "SLC39A10", "IFI27", "ITM2A", "CA4", "CAVIN2", "MT1M", "HERPUD1", "CD320",
              "LGALS3", "SLCO1A2", "SLC7A5", "BSG", "CLDN5", "MT2A", "LY6E", "SLC2A1",
              "NET1", "MT1E", "ADIRF")

features3 = c("MGP", "TIMP1", "CCL2", "SOD2", "SELE", "PDLIM1", "ACKR1", "TM4SF1",
              "NNMT", "EMP1", "VCAM1", "IER3", "ADAMTS9", "IL6", "CD55", "CXCL2",
              "PNP", "ADAMTS1", "DDX21", "YBX3", "PTGS2", "TCIM", "S100A6", "SOCS3",
              "CLU", "SPOCK2", "IFI27", "GNG11", "SLCO2B1", "ADGRF5", "CA4", "SLCO1A2",
              "A2M", "LY6E", "ID1", "ABCG2", "SLC38A5", "SLC39A10", "MFSD2A", "KLF2",
              "CLDN5", "BSG", "HERPUD1", "CD320", "SLC9A3R2", "SLC2A1", "MT1E",
              "CAVIN2", "SLC7A5", "NET1")

features4 = c("SPRY1", "HSPG2", "MGP", "PLVAP", "COL4A1", "SPARC", "COL4A2", "INSR",
              "ANGPT2", "TIMP1", "CALCRL", "CD93", "STC1", "VWA1", "MARCKSL1",
              "IGFBP4", "IGFBP3", "IGFBP7", "HTRA1", "CA2", "ARHGDIB", "MCAM",
              "IVNS1ABP", "VWF", "SOX4", "ABCG2", "IFI6", "SLC38A5", "SLC39A10",
              "TSC22D1", "COBLL1", "SRARP", "CAVIN2", "IFI27", "ITM2A", "CA4",
              "HERPUD1", "CD320", "SLC7A5", "SLCO1A2", "MT1M", "BSG", "LGALS3",
              "CLDN5", "SLC2A1", "LY6E", "NET1", "MT2A", "MT1E", "ADIRF")

Idents(EC) <- "PATH"
plot1 <- DotPlot(EC, features = features1, cols = c("blue", "red"),
                 idents = c("FETALCNS", "TL"), group.by = "PATH",
                 dot.scale = 6, scale = TRUE, scale.by = "size",
                 dot.min = 0, scale.min = 0, scale.max = 100) +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, face = "bold"),
        legend.position = "none") +
  coord_flip() +
  scale_y_discrete(labels = c("FETALCNS" = "Fetal brain", "TL" = "Adult/control brain"))

Idents(EC) <- "PATH2"
plot2 <- DotPlot(EC, features = features2, cols = c("blue", "red"),
                 idents = c("Pathology", "Control"), group.by = "PATH2",
                 dot.scale = 6, scale = TRUE, scale.by = "size",
                 dot.min = 0, scale.min = 0, scale.max = 100) +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, face = "bold"),
        legend.position = "none") +
  coord_flip() +
  scale_y_discrete(labels = c("Pathology" = "Path. brains", 
                              "Control" = "Adult/control brain"))

Idents(EC) <- "PATH3"
plot3 <- DotPlot(EC, features = features3, cols = c("blue", "red"),
                 idents = c("MAL", "CONTROL"), group.by = "PATH3",
                 dot.scale = 6, scale = TRUE, scale.by = "size",
                 dot.min = 0, scale.min = 0, scale.max = 100) +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, face = "bold"),
        legend.position = "none") +
  coord_flip() +
  scale_y_discrete(labels = c("MAL" = "Brain vasc. mal", 
                              "CONTROL" = "Adult/control brain"))

plot4 <- DotPlot(EC, features = features4, cols = c("blue", "red"),
                 idents = c("TUM", "CONTROL"), group.by = "PATH3",
                 dot.scale = 6, scale = TRUE, scale.by = "size",
                 dot.min = 0, scale.min = 0, scale.max = 100) +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, face = "bold")) +
  coord_flip() +
  scale_y_discrete(labels = c("TUM" = "Brain tumours", 
                              "CONTROL" = "Adult/control brain"))

# Combine the plots
combined_plot <- plot1 | plot2 | plot3 | plot4 +plot_layout(guides = "collect")
combined_plot
ggsave("combined_plot.png", plot = combined_plot, width = 20, height = 12)
#------------------------------------------------------------------------------#
