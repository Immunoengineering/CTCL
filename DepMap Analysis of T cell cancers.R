# Molecular Cancer Research Journal: Supporting R code for generating figures
# This script generates plots for analyzing CFL1 and HIF1A CRISPR KO gene effect vs expression
# Corresponding manuscript figures: Figure S3A, S3B and S3C

#Setting work directory
setwd("C:/XXXX") # Adjust to your own path


#Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel) 

#Load gene expression data from DepMap
Expressiondata <- read.csv("OmicsExpression.csv")

#Cell-line metadata from DepMap for a curated list of T cell cancer cell lines
Tcell_cancer_list <- read.csv("T cell cancers.csv")

#=================
#================= Figure S3A: HIF1A CRISPR KO Gene Effect vs Expression =================#

#Loading datasets from DepMap database and renaming the DepMap.ID column as X for ease in merging

#Load CRISPR Data from DepMap for HIF1A
HIF1A_CRISPR <- read.csv("HIF1A CRISPR.csv")
colnames(HIF1A_CRISPR)[colnames(HIF1A_CRISPR) == "Depmap.ID"] <- "X"

#Check HIF1A gene expression vs HIF1A CRISPR KO gene effect across cell lines
HIF1A_data <- inner_join(HIF1A_CRISPR, Expressiondata, by="X")

#renaming columns for 
colnames(HIF1A_data)[2] <- "GeneEffect" 
rownames(HIF1A_data) <- HIF1A_data$X

HIF1A_data$Highlight <- rownames(HIF1A_data) %in% Tcell_cancer_list$DepMap.ID

# Select the Tcell line data
highlighted_points_HIF1A <- HIF1A_data %>%
  filter(Highlight == TRUE)

# Plot for HIF1a gene expression vs HIF1a CRISPR KO gene effect.
ggplot(HIF1A_data, aes(x = GeneEffect, y = HIF1A..3091., color = Highlight)) +
  geom_point(size = 3, aes(alpha = Highlight)) +  
  scale_color_manual(values = c("grey70", "red"), labels = c("Other cell lines", "T cell lines")) +  # Lighter grey for non-highlighted points
  scale_alpha_manual(values = c(0.5, 1)) +  # Setting transparency: 0.5 for non-highlighted, 1 for highlighted
  
  # Using geom_text_repel to avoid overlapping labels
  geom_text_repel(
    data = highlighted_points_HIF1A,  # Using the filtered data for highlighted points
    aes(label = Cell.Line.Name),  # Using Cell.Line.Name column for labels
    size = 3,  # Adjusting label text size
    color = "red",  # Making the label text red
    box.padding = 0.5,  # Padding around the labels
    point.padding = 0.5,  # Padding around the points
    segment.size = 0.2,  # Adjusting the size of the line connecting label to point
    segment.color = "grey50",  # Color of the segment line
    max.overlaps = Inf  # Ensuring that all labels are shown
  ) +
  # Adding a transparent filled orange box highlighting all regions left of -0.85
  annotate("rect", xmin = -Inf, xmax = -0.5, ymin = -Inf, ymax = Inf, 
           fill = "pink", alpha = 0.2) +  # Adjusting alpha for transparency
  
  labs(
    title = "HIF1A CRISPR KO Gene Effect in Cell Lines (DepMap Public Data)",
    x = "Gene Effect by HIF1A CRISPR KO  ",
    y = "HIF1A Expression (Log2TPM+1)",  # Updating y-axis label to match CFL1
    color = "Highlight"  # Legend title
  ) +
  theme_minimal() +  # Using a minimal theme
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))  # Centering the title)

#=================
#================= Figure S3B: CFL1 CRISPR KO Gene Effect vs Expression =================#

#Loading CFL1 datasets from DepMap database and renaming the DepMap.ID column as X for ease in merging


#Loading CRISPR Data from DepMap for CFL1
CFL1_CRISPR <- read.csv("CFL1 CRISPR.csv") 
colnames(CFL1_CRISPR)[colnames(CFL1_CRISPR) == "Depmap.ID"] <- "X"

#CFL1 CRISPR KO Gene effect data combined with gene expression data

CFL1_data <- inner_join(CFL1_CRISPR, Expressiondata, by="X")
colnames(CFL1_data)[2] <- "GeneEffect" 

#Changing rownames to DepMap IDs
rownames(CFL1_data) <- CFL1_data[,1]

#Creating a new column for annotating the cell-line rows as T cell lines

CFL1_data$Highlight <- rownames(CFL1_data) %in% Tcell_cancer_list$DepMap.ID

# Ensure column names are correct and filter for T cell lines
highlighted_points <- CFL1_data %>%
  filter(Highlight == TRUE) %>%
  filter(!is.na(Cell.Line.Name))  # Properly omit NA values in Cell.Line.Name 

# Plot for CFL1 CRISPR Gene effect versus CFL1 gene expression across cell lines
ggplot(CFL1_data, aes(x = GeneEffect, y = CFL1..1072., color = Highlight)) +
  geom_point(size = 3, aes(alpha = Highlight)) +  
  scale_color_manual(values = c("grey70", "red"), labels = c("Other cell lines", "T cell lines")) +  # Lighter grey for non-highlighted points
  scale_alpha_manual(values = c(0.5, 1)) +  # Set transparency: 0.5 for non-highlighted, 1 for highlighted
  
  # Using geom_text_repel to avoid overlapping labels
  geom_text_repel(
    data = highlighted_points,  # Using the filtered data for highlighted points
    aes(label = Cell.Line.Name),  # Using Cell.Line.Name column for labels
    size = 3,  # Adjusting label text size
    color = "red",  # Making the label text red
    box.padding = 0.5,  # Padding around the labels
    point.padding = 0.5,  # Padding around the points
    segment.size = 0.2,  # Adjusting the size of the line connecting label to point
    segment.color = "grey50",  # Color of the segment line
    max.overlaps = Inf  # Ensuring that all labels are shown
  ) +
  # Adding a transparent filled orange box highlighting all regions left of -0.5
  annotate("rect", xmin = -Inf, xmax = -0.5, ymin = -Inf, ymax = Inf, 
           fill = "pink", alpha = 0.2) +  # Adjusting alpha for transparency
  labs(
    title = "CFL1 CRISPR KO Gene Effect in Cell Lines (DepMap Public Data)",
    x = "Gene Effect by CFL1 CRISPR KO",
    y = "CFL1 Expression (Log2TPM+1)",  # Updating y-axis label to match CFL1
    color = "Highlight"  # Legend title
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))  # Centering the title)

#=================
#================= Figure S3C: Density plot for CRISPR KO Gene effect for CFL1 across cell lines =================#


# Creating a new column to distinguish between the two groups (CFL1_data and highlighted_points)
CFL1_data$Group <- "All cell lines"
highlighted_points$Group <- "T cell lines"

# Combining the two datasets for easier plotting
combined_data_CFL1 <- rbind(
  CFL1_data[, c("GeneEffect", "Group")],
  highlighted_points[, c("GeneEffect", "Group")]
)

# Calculating the medians
median_all <- median(combined_data_CFL1$GeneEffect[combined_data_CFL1$Group == "All cell lines"])  # Change "Group1" to your actual group name
median_Tcell <- median(combined_data_CFL1$GeneEffect[combined_data_CFL1$Group == "T cell lines"])  # Change "Group2" to your actual group name


wilcox_test_result_CFL1 <- wilcox.test(CFL1_data$GeneEffect, highlighted_points$GeneEffect)
print(wilcox_test_result_CFL1)

# Calculating Z-score from the p-value
Z <- qnorm(wilcox_test_result_CFL1$p.value / 2) * sign(wilcox_test_result_CFL1$statistic - (length(CFL1_data) * length(highlighted_points)) / 2)

# Calculating the total sample size
N <- length(CFL1_data) + length(highlighted_points)

# Calculating effect size r (Rank-Biserial Correlation)
r <- Z / sqrt(N)

# Creating the combined density plot with a legend
combined_plot_CFL1 <- ggplot(combined_data_CFL1, aes(x = GeneEffect, fill = Group, color = Group)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("blue", "red")) +  # Specify colors for the fill
  scale_color_manual(values = c("blue", "red")) +  # Specify colors for the outline
  labs(
    title = "Combined Density Plot for Gene Effect CFL1 ",
    x = "CRISPR Gene Effect for CFL1 (Chronos)",
    y = "Density",
    fill = "Group",  # Legend title
    color = "Group"  # Legend title for color
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  # Adding vertical dotted lines for medians
  geom_vline(xintercept = median_all, color = "blue", linetype = "dotted", size = 1) +
  geom_vline(xintercept = median_Tcell, color = "red", linetype = "dotted", size = 1) +
  # Adding text for p-value and effect size
  annotate("text", x = min(combined_data_CFL1$GeneEffect), 
           y = max(density(combined_data_CFL1$GeneEffect)$y), 
           label = paste("p =", round(wilcox_test_result_CFL1$p.value, 8), "\nEffect size (r) =", round(r, 3)), 
           hjust = 0, vjust = 1, color = "black", size = 4)

# Displaying the Density plot
print(combined_plot_CFL1)
#==============================================================
