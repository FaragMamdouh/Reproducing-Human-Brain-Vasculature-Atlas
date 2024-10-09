library(Seurat)

# Load the seurat object of the overall merge of all brain sorted CD31+/CD45- endothelial cells 
RPCA <- readRDS(file = "C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
DefaultAssay(RPCA) <- "RNA"
setwd("F:/1- ATLAS/Figure 4")
RPCA@assays[["integrated"]] = NULL
#----------------------------------------------------------------------------#
#Load adult human CNS EC signature genes
AdultCNSsignature <- c("PLP1","CLDN5","GPM6B","TF","BSG","S100B","MBP","GPR85","SLC7A5","PTGDS","SLC2A1","ITM2A","MT3","CNP","SLC38A5","ATP1A2","SPOCK2","SLCO2B1","C1orf61",
                       "APLP1","TUBB4A","QDPR","SCD","MFSD2A","ABCG2","SLC1A3","RN7SK","PPP1R14A","SPOCK3","CRYAB","CKB","TMEM144","PLLP","SRARP","SLC39A10","MAG","ERMN","CLDN11","APOLD1","SCG3",
                       "SCD5","SOX2-OT","TUBA1A","CA4","GPRC5B","CNDP1","CBR1","GPCPD1","LINC00844","CD320")

# Generate mean of expression of the signature genes
AdultCNSsignature_mean <- Matrix::colMeans(RPCA[AdultCNSsignature, ])

# Add matrix to data
RPCA <- AddMetaData(RPCA, metadata = AdultCNSsignature_mean,  col.name = "AdultCNSsignature_mean")
#----------------------------------------------------------------------------#
Idents(RPCA) <- "PATH2"
#to plot the expression of the individual genes of the signature - Figure 4d
pdf("4d.pdf")
DotPlot(RPCA, features = AdultCNSsignature, cols = c("blue", "red"),
        dot.scale = 4, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100)  +coord_flip()
dev.off()

#plot dotplots for expression of the signature in fetal brain vs adult/control
# brain vs pathological brain ECs - Figure 4e (first column)
Idents(RPCA) <- "PATH2"
pdf("4e_first.pdf")
DotPlot(RPCA, features = "AdultCNSsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "CNS signature")
dev.off()

#plot dotplots for expression of the signature by pathology - Figure 4f (first column)
Idents(RPCA) <- "PATH"
pdf("4f_first.pdf")
DotPlot(RPCA, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100)  +
  xlab("")+ ylab("") + scale_x_discrete(labels= "CNS signature")

dev.off()


  #plot dotplots for expression of the signature according to the AV specification - Figure 4g
Idents(RPCA) <- "ECclusters"
levels(RPCA) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_ALL.pdf")
DotPlot(RPCA, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100,) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "ALL")

dev.off()

Idents(RPCA) <- "PATH"

#for fetal brain endothelial cells
FETALCNS <- subset(RPCA, idents = "FETALCNS")
Idents(FETALCNS) <- "ECclusters"
levels(FETALCNS) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_fetal_brain.pdf")
DotPlot(FETALCNS, features = "AdultCNSsignature_mean", 
        cols = c("blue", "red"), dot.scale = 6, 
        scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100)  +
  xlab("")+ ylab("") + scale_x_discrete(labels= "Fetal brain")
dev.off()

remove(FETALCNS)

#for adult/control brain endothelial cells
TL <- subset(RPCA, idents = c("TL"))
Idents(TL) <- "ECclusters"
levels(TL) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_adult_control_brain.pdf")
DotPlot(TL, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100)  +
  xlab("")+ ylab("") + scale_x_discrete(labels= "Adult/control brain")
dev.off()
remove(TL)
#for brain vascular malformation (arteriovenous malformation) endothelial cells (AVM/Brain vasc. MAL)
AVM <- subset(RPCA, idents = c("AVM"))
Idents(AVM) <- "ECclusters"
levels(AVM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_Brain vasc. mal.pdf")
DotPlot(AVM, features = "AdultCNSsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "Brain vasc. mal")
dev.off()
remove(AVM)

#for lower-grade glioma (LGG) endothelial cells
LGG <- subset(RPCA, idents = c("LGG"))
Idents(LGG) <- "ECclusters"
levels(LGG) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_LGG.pdf")
DotPlot(LGG, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, scale.min = 0, 
        scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "LGG")
dev.off()
remove(LGG)
#for high-grade glioma/glioblastoma (GBM) endothelial cells
GBM <- subset(RPCA, idents = c("GBM"))
Idents(GBM) <- "ECclusters"
levels(GBM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_GBM.pdf")
DotPlot(GBM, features = "AdultCNSsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "GBM")
dev.off()
remove(GBM)
#for metastasis (MET) endothelial cells
MET <- subset(RPCA, idents = c("MET"))
Idents(MET) <- "ECclusters"
levels(MET) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_MET.pdf")

DotPlot(MET, features = "AdultCNSsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "MET")
dev.off()
remove(MET)
#for meningioma (MEN) endothelial cells
MEN <- subset(RPCA, idents = c("MEN"))
Idents(MEN) <- "ECclusters"
levels(MEN) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_MEN.pdf")
DotPlot(MEN, features = "AdultCNSsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "MEN")
dev.off()
remove(MEN)

#for all pathological brain endothelial cells (PATH)
Idents(RPCA) <- "PATH2"
PATH <- subset(RPCA, idents = c("Pathology"))
Idents(PATH) <- "ECclusters"
levels(PATH) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_PATH.pdf")
DotPlot(PATH, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "PATH")
dev.off()
remove(PATH)


#for all brain tumor endothelial cells (Brain tumors)
Idents(RPCA) <- "PATH3"
TUM <- subset(RPCA, idents = c("TUM"))
Idents(TUM) <- "ECclusters"
levels(TUM) <- c("Proliferating Stem-to-EC","Stem-to-EC", "Proliferating cell", "Proliferating EndoMT","EndoMT","Mitochondrial", "Large vein", "Vein", "Venule", "Capillary","Angiogenic capillary", "Arteriole", "Artery","Large artery")
pdf("4g_Brain_tum.pdf")
DotPlot(TUM, features = "AdultCNSsignature_mean", cols = c("blue", "red"), 
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100) +
  xlab("")+ ylab("") + scale_x_discrete(labels= "Brain tumours")
dev.off()
remove(TUM)

Idents(RPCA) <- "PATH"
#------------------------------------------------------------------------------# 
#violn plot for expression of some CNS specific genes - Figure 4h
Idents(RPCA) <- "PATH2"
levels(RPCA) <- c("Fetal","Control","Pathology")

library(gridExtra)
genes <- c("SPOCK3", "GPCPD1", "BSG", "PPP1R14A", "CD320", "SLC38A5", 
           "ADGRA2", "SLC2A1", "CLDN5", "TJP1", "OCLN")

plots <- list()
for (gene in genes) {
  p <- VlnPlot(RPCA, features = c(gene), log = TRUE, pt.size = 0, 
               cols = c('Fetal' = '#A0E4A8', 'Control' = 'turquoise3',
                        'Pathology' = '#f8766d'))+ NoLegend()+ ylab("") +xlab("")
  plots[[gene]] <- p
}
pdf("4h.pdf")
grid.arrange(grobs = plots, ncol = 3)

dev.off()
#------------------------------------------------------------------------------#

#Load adult human peripheral EC signature genes
AdultPeriphsignature <-
  c("LDHA","NCOA7","TM4SF1","NNMT","CXCL3","FABP5","SLPI","GAS5","IGFBP4","DDX21","POSTN","IL33","MT-RNR2","CD36","THBS1","PDLIM1","HIF1A","RND1","CD93","SLCO2A1",
    "ICAM1","CRHBP","TFPI","CALCRL","COL15A1","SOCS3","PNP","SOD2","ADAMTS9","HLA-DRB5","CCL2","CXCL2","FABP4","PLVAP","PLPP3","SERPINE1","ACKR1","DNASE1L3","FCN3","THBD",
    "IGKC","STC1","IGFBP5","HMOX1","EMP1","AQP1","CSF3","IL6","CXCL8","SELE")

# Generate mean of expression of the signature genes
AdultPeriphsignature_mean <- Matrix::colMeans(RPCA[AdultPeriphsignature, ])

# Add matrix to data
RPCA <- AddMetaData(RPCA, metadata = AdultPeriphsignature_mean,  col.name = "AdultPeriphsignature_mean")

#plot dotplots for expression of the signature by pathology - Figure 4f (second column)
pdf("4f_second.pdf")
Idents(RPCA) <- "PATH"
DotPlot(RPCA, features = "AdultPeriphsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0,
        scale.min = 0, scale.max = 100)  +
  xlab("") + ylab("") + scale_x_discrete(labels = "Peripheral Signature")

dev.off()

#plot dotplots for expression of the signature in fetal brain vs adult/control brain vs pathological brain ECs - Figure 4e (second column)
pdf("4e_second.pdf")
Idents(RPCA) <- "PATH2"
DotPlot(RPCA, features = "AdultPeriphsignature_mean", cols = c("blue", "red"),
        dot.scale = 6, scale = TRUE, scale.by = "size", dot.min = 0, 
        scale.min = 0, scale.max = 100) +
  xlab("") + ylab("") + scale_x_discrete(labels = "Peripheral Signature")
dev.off()
Idents(RPCA) <- "ECclusters"
ggplot = FeaturePlot(RPCA, 
            features = c("AdultCNSsignature_mean","AdultPeriphsignature_mean"),
            split.by = "PATH2", blend = TRUE, reduction = 'UMAP', 
            label = TRUE,
            repel = TRUE, label.size = 4,pt.size = 0.1,
            min.cutoff = "q05", max.cutoff = "q90",raster=FALSE,
            blend.threshold = 0.5, combine = FALSE)

Pathological <- ggplot[[3]] + ggtitle("Pathological brains") +
  theme(legend.position = "none")

Adult <-ggplot[[7]] + ggtitle("Adult/control brains") +
  theme(legend.position = "none")
FETAL <- ggplot[[11]] +ggtitle("Fetal brains") +
  theme(legend.position = "none")
ggsave("Pathological_brains.pdf", plot = Pathological, width = 7, height = 5, units = "in")

ggsave("Adult_control_brains.pdf", plot = Adult, width = 7, height = 5, units = "in")

ggsave("Fetal_brains.pdf", plot = FETAL, width = 7, height = 5, units = "in")
