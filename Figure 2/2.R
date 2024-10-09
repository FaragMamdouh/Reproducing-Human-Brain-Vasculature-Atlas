
library(readxl)
library(VennDiagram)
library(ggplot2)
library(dplyr)
all <- "F:/1- ATLAS/supp tables/Supplementary Table 17.xlsx"
setwd("F:/1- ATLAS/Figure 2/mouse and human")
#read supplementary Table 17 file
sheet_names <- excel_sheets(all)

# Retrieve the tissues from the Excel file

for (sheet in sheet_names) {
  data <- read_excel(all, sheet = sheet)
  venn.plot <- venn.diagram(
    x = list(Set1 = na.omit(data$human), Set2 = na.omit(data$mouse)),  # Adjust based on your column names
    filename = NULL,
    fill = c("#c6cdea", "#f1c6c7"),     
    main = paste( sheet),  
    col = c("#7390cd", "#e79eb0"),
    main.fontface = "bold",  
    main.cex = 1.5,  
    cat.cex = 1.2, 
    cex = 1.2, 
    lwd = 2  ,
    category.names = c("Human", "Mouse")
  )
  
  grid::grid.draw(venn.plot)
  dev.off()
  ggsave(filename = paste0(gsub(" ", "_", sheet), ".tiff"),
         plot = venn.plot,  width = 10, height = 8, units = "in")
}
# Create a data frame summarizing common and specific counts for different clusters based on venn diagram

bottom <- data.frame(
  Cluster = c("Large artery", "Artery", "Arteriole", "Capillary", "Venule", "Large vein"),
  Common = c(87, 32, 30, 2, 1, 47),
  Human_specific = c(181, 72, 248, 67, 75, 276),
  Mouse_specific = c(299, 265, 79, 11, 64, 153)
)

# Calculate the total counts for each cluster

bottom$total <- rowSums(bottom[,-1])

# Calculate percentages for common, human-specific, and mouse-specific counts
bottom$Common_percent <- (bottom$Common / bottom$total) * 100
bottom$Human_specific_percent <- (bottom$Human_specific / bottom$total) * 100
bottom$Mouse_specific_percent <- (bottom$Mouse_specific / bottom$total) * 100

# Transform the data frame to a long format for plotting
bottom <- bottom %>%
  select(Cluster, Common_percent, Human_specific_percent, Mouse_specific_percent) %>%
  pivot_longer(cols = -Cluster, names_to = "Type", values_to = "Percentage")

# Clean up the 'Type' column for better readability
bottom$Type <- gsub("_percent", "", bottom$Type)
bottom$Type <- gsub("_", " ", bottom$Type)

# Create a stacked bar plot of the percentages
plot <- ggplot(bottom, aes(x = Cluster, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(y = "Percentage", fill = "", x = "Brain EC clusters") +
  scale_fill_manual(values = c("Common" = "lightgrey", "Human specific" = "lightblue", "Mouse specific" = "lightcoral")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 100)) 

ggsave("cluster_percentages.png", plot = plot, width = 8, height = 5, dpi = 300)

