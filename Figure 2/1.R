library("heatmaply")
library(Seurat)
library(dplyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
setwd("C:/Users/Farag/Downloads/GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")
list.files()
EC <- readRDS(file = "GSE256493_Overall_merge_of_all_brain_sorted_endothelial_cells_seurat_object.rds")

setwd("F:/1- ATLAS/Figure 2")
png("2a1.png")
DimPlot(EC, reduction = "UMAP", label = TRUE, repel = TRUE, label.size = 4,
        pt.size = 1, group.by = "ECclusters", raster=FALSE,
        cols = c('Large artery' = '#a1000b','Artery' = '#d0000d', 
                 'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe', 
                 'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd', 
                 'Vein' = '#29e2e6','Large vein' = '#0000b3', 
                 'Mitochondrial' = '#f1b77c', 'EndoMT' = '#70fb52', 
                 'Proliferating EndoMT' = '#affd2d', 'Proliferating cell' = '#FBFF00', 
                 'Stem-to-EC' = '#4cd7a4', 'Proliferating Stem-to-EC' = '#89e56b')) +
  labs(title = "Overall merge brain ECs") +
  theme(title = element_text(size = 12, face = "plain", hjust = .5))
dev.off()

Idents(EC) <- "PATH2"
EC$PATH2 <- factor(EC$PATH2, levels = c("Fetal", "Control", "Pathology"))

# Figure 2a2
png("2a2.png")
DimPlot(EC, reduction = "UMAP", pt.size = 1, 
        label = FALSE, repel = TRUE, split.by = "PATH2", label.size = 4, raster = FALSE, 
        cols = c("green", "cyan", "red")) + 
  NoLegend() + 
  facet_wrap(~PATH2, labeller = as_labeller(c(
    "Fetal" = "Fetal Brain ECs",
    "Control" = "Adult/Control Brain ECs",
    "Pathology" = "Pathological Brain ECs"
  )))  
dev.off()
#------------------------------Figure 2B---------------------------------------#
EC <- AddMetaData(object = EC, metadata = "EC" , col.name = "EC")
Idents(EC) <- "ECclusters"
levels(EC) <- c("Large artery", "Artery", "Arteriole","Angiogenic capillary", "Capillary", "Venule", "Vein", "Large vein","Mitochondrial", "EndoMT","Proliferating EndoMT", "Proliferating cell","Stem-to-EC","Proliferating Stem-to-EC")

#extract metadata to do compostional bargraphs
metadata <- as.data.frame(as.matrix(EC@meta.data))

#arrange metadadata of EC AV clusters
metadata$ECclusters <- factor(metadata$ECclusters, 
                              levels=c('Large artery', 'Artery', 'Arteriole',
                                       'Angiogenic capillary', 'Capillary',
                                       'Venule', 'Vein', 'Large vein',
                                       'Mitochondrial', 'EndoMT', 
                                       'Proliferating EndoMT', 'Proliferating cell',
                                       'Stem-to-EC', 'Proliferating Stem-to-EC'))

#arrange metadadata of EC AV sub-clusters
metadata$ECsubclusters <- factor(metadata$ECsubclusters, 
                                 levels=c('Large artery','Artery1','Artery2',
                                          'Arteriole1','Arteriole2','Capillary1',
                                          'Capillary2','Angiogenic capillary1',
                                          'Angiogenic capillary2','Venule1',
                                          'Venule2','Vein1','Vein2','Large vein',
                                          'EndoMT1','EndoMT2','Proliferating EndoMT',
                                          'Proliferating cell1','Proliferating cell2',
                                          'Proliferating Stem-to-EC','Stem-to-EC1',
                                          'Stem-to-EC2','Mitochondrial'))

#arrange meradadata of Fetal vs Control vs Pathology annotation
metadata$PATH2 <- factor(metadata$PATH2, levels=c('Fetal','Control','Pathology'))


#arrange meradadata of Patients
metadata$Patient <- factor(metadata$Patient,levels=c('FETALCNS1','FETALCNS2-3','FETALCNS4','FETALCNS5','TL1','TL2','TL3','TL4','TL5','TL6','TL7','TL9','TL10','AVM1','AVM2','AVM3','AVM4','AVM5','LGG1','LGG2','LGG3','LGG4','LGG5','LGG6','GBM1','GBM2','GBM3','GBM4','GBM5','GBM6','GBM7','GBM8','MET1','MET2','MET3','MET4','MET5',
                                                     'MEN1','MEN2','MEN3','MEN4','MEN5'))

#----------------------------Figure 2b-d---------------------------------------#
cols = c('Large artery' = '#a1000b','Artery' = '#d0000d', 
         'Arteriole' = '#d13b18', 'Angiogenic capillary' = '#c300fe',
         'Capillary' = '#fb00c9', 'Venule' = '#6dc4fd', 'Vein' = '#29e2e6',
         'Large vein' = '#0000b3', 'Mitochondrial' = '#f1b77c', 
         'EndoMT' = '#70fb52', 'Proliferating EndoMT' = '#affd2d', 
         'Proliferating cell' = '#FBFF00', 'Stem-to-EC' = '#4cd7a4',
         'Proliferating Stem-to-EC' = '#89e56b')


#----------------------------ALL cells
ALLCELLS <- ggplot(metadata, aes(x = metadata[["EC"]], fill = metadata[["ECclusters"]])) + 
  geom_bar(position = "fill") + 
  labs(fill = NULL) + 
  xlab("All Cells") + 
  theme_void() + 
  scale_fill_manual(values = cols) + 
  theme(axis.text.y = element_text(face = "bold"), legend.position = "none") + 
  scale_x_discrete(labels = c("EC" = "All cells")) + 
  coord_flip() + 
  scale_y_reverse() +
  scale_y_continuous(limits = c(0, 1))
#---------------------------Fetal Brain and Control
Fetal_control <- metadata %>% filter(PATH2 %in% c("Fetal", "Control"))
Fetal_control_plot <- ggplot(Fetal_control, aes(x= Fetal_control[["PATH2"]], 
                                                fill= Fetal_control[["ECclusters"]])) + 
  geom_bar(position = "fill")+ 
  labs(fill= NULL)+ 
  xlab("Cluster")+ 
  ylab("") + 
  scale_fill_manual(values=cols) +
  coord_flip() +
  theme_void() +
  theme( axis.text.y = element_text(face = "bold"), legend.position = "none")+
  scale_x_discrete(labels = c("Fetal" = "Fetal Brain", 
                              "Control" = "TL"))+
  scale_y_reverse()+
  scale_y_continuous(limits = c(0, 1))


#---------------------------patholical
patholical <- metadata %>% filter(PATH2 %in% c("Pathology"))
patholicalplot <-  ggplot(patholical, aes(x= patholical[["PATH2"]], 
                           fill= patholical[["ECclusters"]])) + 
  geom_bar(position = "fill")+ 
  labs(fill= NULL)+ 
  xlab("Cluster")+ 
  ylab("") + 
  scale_fill_manual(values=cols) +
  coord_flip() +
  theme_void() +
  theme( axis.text.y = element_text(face= "bold"))+
  scale_x_discrete(labels = c("Pathology" = "All pathological brains"))+
  scale_y_reverse()+
  scale_y_continuous(limits = c(0, 1))


#---------------------------AVM
AVM <- metadata%>%
  filter(PATH %in% "AVM")
AVMplot =  ggplot(AVM, aes(x= AVM[["PATH"]], fill= AVM[["ECclusters"]])) + 
  geom_bar(position = "fill")+ 
  labs(fill= NULL)+ xlab("Cluster")+ 
  ylab("") + 
  scale_fill_manual(values=cols)+
  coord_flip() +
  theme_void() +
  theme( axis.text.y = element_text(face = "bold")
         , legend.position = "none")+
  scale_x_discrete(labels = c("AVM" = "Brain vascular malformattion" ))+
  scale_y_reverse() +
  scale_y_continuous(limits = c(0, 1))

#---------------------------Brain Tumors
BrainTum <- metadata%>%
  filter(PATH3 %in% "TUM")
BrainTumplot =  ggplot(BrainTum, aes(x= BrainTum[["PATH3"]], fill= BrainTum[["ECclusters"]])) + 
  geom_bar(position = "fill")+ 
  labs(fill= NULL)+ xlab("Cluster")+ 
  ylab("") + 
  scale_fill_manual(values=cols)+
  coord_flip() +
  theme_void() +
  theme( axis.text.y = element_text(face = "bold"), legend.position = "none")+
  scale_x_discrete(labels = c("TUM" = "Brain tumors" ))+
  scale_y_reverse() +
  scale_y_continuous(limits = c(0, 1))

#---------------------------TUMOrs_separate
Tumors_separate <- metadata %>% filter(PATH %in% c("GBM", "LGG", "MEN", "MET"))

Tumors_separate_plot <- ggplot(Tumors_separate, aes(x= PATH, fill= ECclusters)) + 
  geom_bar(position = "fill") + 
  labs(fill= NULL) + 
  xlab("Cluster") + 
  ylab("") + 
  scale_fill_manual(values=cols) + 
  coord_flip() +
  theme_void() +
  theme(axis.text.y = element_text(face= "bold"), legend.position = "none") +
  scale_y_reverse() +
  scale_y_continuous(limits = c(0, 1))



library("patchwork")
combined_percentage_2 <- (ALLCELLS / Fetal_control_plot / patholicalplot / 
                            AVMplot / BrainTumplot / Tumors_separate_plot)
ggsave("combined_percentage2.png", 
       combined_percentage_2,units = "cm", height = 15)
#------------------------------------------------------------------------------#


# Input the the top AV specification features to plot
features = c("MGP","FBLN5","LTBP4","S100A6","ELN",
             "GLUL","ADAMTS1","SEMA3G","ALPL","GJA4",
             "SLC2A1","LGALS3","CD320","AIF1L","NET1",
             "ESM1","ANGPT2","APLN","PLVAP","CA2",
             "CA4","SLCO1A2","BSG","MFSD2A","SLC16A1",
             "IL1R1","PRCP","RAMP3","PRSS23","HYAL2",
             "POSTN","PTGDS","DNASE1L3","NNMT","EFEMP1",
             "SELE","ACKR1","IL6","CCL2","VCAM1",
             "APOE","TAGLN","ACTA2","NDUFA4L2","RGS5",
             "THY1","TPM2","LGALS1","COL1A2","COL1A1",
             "TTN","NPIPB5","SMG1","MT-CO3","MT-RNR1",
             "CENPF","TOP2A","STMN1","MKI67","PTTG1",
             "HIST1H4C","HMGN2","BEX1","HMGB2","UBE2C",
             "MT3","FABP7","VEGFA","PTPRZ1","SEC61G")

# Get the average expression per cluster and convert to data frame
AV_markers <- AverageExpression(EC, assays = "RNA", features = features, 
                                return.seurat = FALSE, group.by = "ECclusters",
                                add.ident = NULL, slot = "data", verbose = TRUE)


markers_dataframe <- as.data.frame(AV_markers)


# Arrange the dataframe according to the order you need
AV_markers_dataframe2 <- markers_dataframe[, c("RNA.Large.artery","RNA.Artery", "RNA.Arteriole","RNA.Angiogenic.capillary", "RNA.Capillary","RNA.Venule","RNA.Vein", "RNA.Large.vein", "RNA.EndoMT","RNA.Proliferating.EndoMT","RNA.Mitochondrial", "RNA.Proliferating.cell","RNA.Proliferating.Stem.to.EC","RNA.Stem.to.EC")]


heatmaply(AV_markers_dataframe2, scale = "row",
          colors = c("#0571b0","white","#ca0020"), 
          grid_gap = 1,dendrogram = ("none")) %>% 
  layout(width=1000, height = 1400)

