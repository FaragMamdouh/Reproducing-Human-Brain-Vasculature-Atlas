library("Seurat")
setwd("F:/1- ATLAS/Figure 5")

#Load MHC class II signature genes
MHCCLASSII <- c("CD74","HLA-DRB5","HLA-DRB1","HLA-DQB2","HLA-DQB1","HLA-DPB1",
                "HLA-DOB","HLA-DMB","HLA-DRA","HLA-DQA2","HLA-DQA1","HLA-DPA1",
                "HLA-DOA","HLA-DMA")  

# Generate mean of expression of the signature genes
MHCCLASSII_mean <- Matrix::colMeans(RPCA[MHCCLASSII, ])

# Add matrix to data
RPCA <- AddMetaData(RPCA, metadata = MHCCLASSII_mean, col.name = "MHCCLASSII_mean")
Idents(RPCA) <- "PATH2"
#to plot the expression of the individual genes of the signature - Figure 5d
pdf("5d_first.pdf")
DotPlot(RPCA, features = MHCCLASSII, cols = c("blue", "red"), dot.scale = 6,
        scale = TRUE, scale.by = "size", dot.min = 0, scale.min = 0, 
        scale.max = 100, col.min = 0, col.max = 10)+coord_flip() +
  xlab("") + ylab("") 
dev.off()
#plot dotplots for expression of the signature in fetal brain vs adult/control 
# brain vs pathological brain ECs - Figure 5e (first column)
Idents(RPCA) <- "PATH2"
pdf("5e_first.pdf")
DotPlot(RPCA, features = "MHCCLASSII_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")
dev.off()
#plot dotplots for expression of the signature by pathology - Figure 5f
Idents(RPCA) <- "PATH"
pdf("5f.pdf")
DotPlot(RPCA, features = "MHCCLASSII_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()


#plot dotplots for expression of the signature according to the AV specification 
# - Figure 5g

#for all brain endothelial cells (ALL cells)
Idents(RPCA) <- "ECclusters"
levels(RPCA) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                  "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                  "Vein", "Venule", "Capillary","Angiogenic capillary", 
                  "Arteriole", "Artery","Large artery")
pdf("5g_all.pdf")
DotPlot(RPCA, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()

Idents(RPCA) <- "PATH"
#for fetal brain endothelial cells
FETALCNS <- subset(RPCA, idents = c("FETALCNS"))
Idents(FETALCNS) <- "ECclusters"
levels(FETALCNS) <- c("Proliferating Stem-to-EC","Stem-to-EC", 
                      "Proliferating cell", "Proliferating EndoMT","EndoMT",
                      "Mitochondrial", "Large vein", "Vein", "Venule",
                      "Capillary","Angiogenic capillary", "Arteriole",
                      "Artery","Large artery")
pdf("5g_fetal.pdf")
DotPlot(FETALCNS, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", 
        dot.min = 0, scale.min = 0, scale.max = 100, 
        col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")
dev.off()
remove(FETALCNS)
#for adult/control brain endothelial cells
TL <- subset(RPCA, idents = c("TL"))
Idents(TL) <- "ECclusters"
levels(TL) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", 
                "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                "Vein", "Venule","Capillary","Angiogenic capillary","Arteriole", 
                "Artery","Large artery")

pdf("5g_TL.pdf")
DotPlot(TL, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(TL)
#for brain vascular malformation (arteriovenous malformation) endothelial cells 
# (AVM/Brain vasc. MAL)
AVM <- subset(RPCA, idents = c("AVM"))
Idents(AVM) <- "ECclusters"
levels(AVM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", 
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", 
                 "Vein", "Venule", "Capillary","Angiogenic capillary", 
                 "Arteriole", "Artery","Large artery")
pdf("5g_MAL.pdf")
DotPlot(AVM, features = "MHCCLASSII_mean",cols = c("blue", "red"),dot.scale = 6,
        scale = TRUE, scale.by = "size", dot.min = 0, scale.min = 0, 
        scale.max = 100, col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(AVM)
#for lower-grade glioma (LGG) endothelial cells
LGG <- subset(RPCA, idents = c("LGG"))
Idents(LGG) <- "ECclusters"
levels(LGG) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                 "Vein", "Venule", "Capillary","Angiogenic capillary", 
                 "Arteriole", "Artery","Large artery")
pdf("5g_LGG.pdf")
DotPlot(LGG, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10) +
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")
dev.off()
remove(LGG)
#for high-grade glioma/glioblastoma (GBM) endothelial cells
GBM <- subset(RPCA, idents = c("GBM"))
Idents(GBM) <- "ECclusters"
levels(GBM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", 
                 "Vein", "Venule", "Capillary","Angiogenic capillary", 
                 "Arteriole", "Artery","Large artery")
pdf("5g_GBM.pdf")
DotPlot(GBM, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10)+
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(GBM)
#for metastasis (MET) endothelial cells
MET <- subset(RPCA, idents = c("MET"))
Idents(MET) <- "ECclusters"
levels(MET) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                 "Vein", "Venule", "Capillary","Angiogenic capillary", 
                 "Arteriole", "Artery","Large artery")
pdf("5g_MET.pdf")

DotPlot(MET, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10)+
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(MET)
#for meningioma (MEN) endothelial cells
MEN <- subset(RPCA, idents = c("MEN"))
Idents(MEN) <- "ECclusters"
levels(MEN) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", 
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                 "Vein", "Venule", "Capillary","Angiogenic capillary",
                 "Arteriole", "Artery","Large artery")
pdf("5g_MEN.pdf")
DotPlot(MEN, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10)+
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(MEN)
#for all pathological brain endothelial cells (PATH)
Idents(RPCA) <- "PATH2"
PATH <- subset(RPCA, idents = c("Pathology"))
Idents(PATH) <- "ECclusters"
levels(PATH) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                  "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                  "Vein", "Venule", "Capillary","Angiogenic capillary", 
                  "Arteriole", "Artery","Large artery")
pdf("5g_PATH.pdf")

DotPlot(PATH, features = "MHCCLASSII_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100, col.min = 0, col.max = 10)+
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()
remove(PATH)
#for all brain tumor endothelial cells (Brain tumors)
Idents(RPCA) <- "PATH3"
TUM <- subset(RPCA, idents = c("TUM"))
Idents(TUM) <- "ECclusters"
levels(TUM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell",
                 "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein",
                 "Vein", "Venule", "Capillary","Angiogenic capillary", 
                 "Arteriole", "Artery","Large artery")
pdf("5g_brain_tum.pdf")

DotPlot(TUM, features = "MHCCLASSII_mean",cols = c("blue", "red"),dot.scale = 6,
        scale = TRUE, scale.by = "size", dot.min = 0, scale.min = 0,
        scale.max = 100, col.min = 0, col.max = 10)+
  xlab("")+ylab("")+  scale_x_discrete(labels = "MHC class II signature")

dev.off()

#to plot a UMAP with the blend of the MHC class II and Adult CNS signature in fetal brain vs 
# adult/control brain vs pathological brain ECs - Figure 5a-c
Idents(RPCA) <- "ECclusters"
ggplot <- FeaturePlot(RPCA, features = c("MHCCLASSII_mean","AdultCNSsignature_mean"),
            split.by = "PATH2", blend = TRUE, reduction = 'UMAP',label = TRUE, 
            repel = TRUE, label.size = 4,pt.size = 0.1, 
            min.cutoff = "q05", max.cutoff = "q90",
            raster=FALSE,blend.threshold = 0.5, combine = FALSE)
Pathological <- ggplot[[3]] + ggtitle("Pathological brains") +
  theme(legend.position = "none")

Adult <-ggplot[[7]] + ggtitle("Adult/control brains") +
  theme(legend.position = "none")
FETAL <- ggplot[[11]] +ggtitle("Fetal brains") +
  theme(legend.position = "none")
ggsave("Pathological_brains.pdf", plot = Pathological, width = 7, height = 5, units = "in")

ggsave("Adult_control_brains.pdf", plot = Adult, width = 7, height = 5, units = "in")

ggsave("Fetal_brains.pdf", plot = FETAL, width = 7, height = 5, units = "in")
