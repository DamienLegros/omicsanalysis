### Author: Damien Legros
### Date: 05/12/2024
### Comment: Small Data Analysis for Integrative Microbiome and Metabolomics

###############################################
### Add packages necessary for the analysis ###
###############################################

library(phyloseq)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ape)

####################
### Data loading ###
####################

# Set working directory
setwd("omicsanalysis")

# Load the .RData file
load("metagenomics10y.RData")

# To make the iTOL tree
tree <- phy_tree(phyX)
write.tree(tree, file = "phylogenetic_tree.nwk")

#################
### Functions ###
#################

# Function to check if a row has constant values
is_constant <- function(x) {
  all(x == x[1], na.rm = TRUE)  # Check if all values in the row are the same
}

###########################
### Data pre-processing ###
###########################

### Sample data preprocessing

# Create the dataframe
df_sample_data <- data.frame(sample_data(phyX)) 

# Take out the Timepoint column
df_sample_data <- df_sample_data %>% select(-Timepoint)

# Remove the NA samples
df_sample_data <- na.omit(df_sample_data)

# Change the type to num
df_sample_data <- df_sample_data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

## Select overall variables
df_sample_data_overall <- df_sample_data[, 1:2]

## Select food variables
df_sample_data_food <- df_sample_data[, 3:42]

# Take out the abundance values that are constant
df_sample_data_food <- df_sample_data_food[, !apply(df_sample_data_food, 2, is_constant)]

## Select nutrients variables
df_sample_data_nutrients <- df_sample_data[, 43:158]

# Take out the abundance values that are constant
df_sample_data_nutrients <- df_sample_data_nutrients[, !apply(df_sample_data_nutrients, 2, is_constant)]

## Select metabolomics variables
df_sample_data_metabolomics <- df_sample_data[, 159:239]

# Take out the abundance values that are constant
df_sample_data_metabolomics <- df_sample_data_metabolomics[, !apply(df_sample_data_metabolomics, 2, is_constant)]

### OTU table preprocessing

# Convert the OTU table to dataframe 
df_otu_table <- data.frame(otu_table(phyX))

# Normalize OTU table by column sums
df_otu_table <- df_otu_table %>% mutate(across(everything(), ~ . / sum(.)))

# Convert the taxas table to dataframe 
df_tax_table <- data.frame(tax_table(phyX))

## Group by Phylum and sum OTU counts
df_otu_table$Phylum <- df_tax_table$Phylum
df_otu_table_phylum <- df_otu_table %>% group_by(Phylum) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
df_otu_table_phylum <- data.frame(df_otu_table_phylum, row.names = 1)
df_otu_table_phylum$Phylum <- NULL
df_otu_table$Phylum <- NULL

# Take out the abundance values that are constant across all samples and transpose
df_otu_table_phylum <- df_otu_table_phylum[!apply(df_otu_table_phylum, 1, is_constant), ]
df_otu_table_phylum <- t(df_otu_table_phylum)

# Keep only the sample rows with data
df_otu_table_phylum <- df_otu_table_phylum[rownames(df_sample_data), ]

## Group by Class and sum OTU counts
df_otu_table$Class <- df_tax_table$Class
df_otu_table_class<- df_otu_table %>% group_by(Class) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
df_otu_table_class <- data.frame(df_otu_table_class, row.names = 1)
df_otu_table_class$Class <- NULL
df_otu_table$Class <- NULL

# Take out the abundance values that are constant across all samples and transpose
df_otu_table_class <- df_otu_table_class[!apply(df_otu_table_class, 1, is_constant), ]
df_otu_table_class <- t(df_otu_table_class)

# Keep only the sample rows with data
df_otu_table_class <- df_otu_table_class[rownames(df_sample_data), ]

## Group by Order and sum OTU counts
df_otu_table$Order <- df_tax_table$Order
df_otu_table_order<- df_otu_table %>% group_by(Order) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
df_otu_table_order <- data.frame(df_otu_table_order, row.names = 1)
df_otu_table_order$Order <- NULL
df_otu_table$Order <- NULL

# Take out the abundance values that are constant across all samples and transpose
df_otu_table_order <- df_otu_table_order[!apply(df_otu_table_order, 1, is_constant), ]
df_otu_table_order <- t(df_otu_table_order)

# Keep only the sample rows with data
df_otu_table_order <- df_otu_table_order[rownames(df_sample_data), ]

## Group by Family and sum OTU counts
df_otu_table$Family <- df_tax_table$Family
df_otu_table_family <- df_otu_table %>% group_by(Family) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
df_otu_table_family <- data.frame(df_otu_table_family, row.names = 1)
df_otu_table_family$Family <- NULL
df_otu_table$Family <- NULL

# Take out the abundance values that are constant across all samples and transpose
df_otu_table_family <- df_otu_table_family[!apply(df_otu_table_family, 1, is_constant), ]
df_otu_table_family <- t(df_otu_table_family)

# Keep only the sample rows with data
df_otu_table_family <- df_otu_table_family[rownames(df_sample_data), ]

#####################
### Data analysis ###
#####################

### Correlation matrices

# Open the device to show the plots
dev.new()

# Use the default color palette from pheatmap
default_palette <- colorRampPalette(c("#5f84bc", "#fdfec0", "#d6402e"))(100)

## Only sample data correlation matrices

# Overall-Food Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_overall), as.matrix(df_sample_data_food), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), main = "Overall-Food Correlations", filename = "plots/overall_food_corr.png")

# Overall-Nutrients Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_overall), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), fontsize_col = 6, main = "Overall-Nutrients Correlations", filename = "plots/overall_nutrients_corr.png")

# Overall-Metabolomics Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_overall), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), fontsize_col = 6, main = "Overall-Metabolomics Correlations", filename = "plots/overall_metabolomics_corr.png")

# Food-Nutrients Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_food), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), fontsize_col = 6, main = "Food-Nutrients Correlations", filename = "plots/food_nutrients_corr.png")

# Food-Metabolomics Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_food), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), fontsize_col = 6, main = "Food-Metabolomics Correlations", filename = "plots/food_metabolomics_corr.png")

# Nutrients-Metabolomics Correlation matrix heatmap

cor_matrix <- cor(as.matrix(df_sample_data_nutrients), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-1, 1, length.out = length(default_palette) + 1), fontsize_col = 6, fontsize_row = 6, main = "Nutrients-Metabolomics", filename = "plots/nutrients_overall_corr.png")

## Phylum Correlation matrices

# Phylum-Overall Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_phylum), as.matrix(df_sample_data_overall), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), main = "Phylum-Overall Correlations", filename = "plots/phylum_overall_corr.png")

# Phylum-Food Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_phylum), as.matrix(df_sample_data_food), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), main = "Phylum-Food Correlations", filename = "plots/phylum_food_corr.png")

# Phylum-Nutrients Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_phylum), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_col = 4, main = "Phylum-Nutrients Correlations", filename = "plots/phylum_nutrients_corr.png")

# Phylum-Metabolomics Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_phylum), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_col = 5, main = "Phylum-Metabolomics Correlations", filename = "plots/phylum_metabolomics_corr.png")

## Class Correlation matrices

# Class-Overall Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_class), as.matrix(df_sample_data_overall), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), main = "Class-Overall Correlations", filename = "plots/class_overall_corr.png")

# Class-Food Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_class), as.matrix(df_sample_data_food), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), main = "Class-Food Correlations", filename = "plots/class_food_corr.png")

# Class-Nutrients Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_class), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_col = 4, main = "Class-Nutrients Correlations", filename = "plots/class_nutrients_corr.png")

# Class-Metabolomics Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_class), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_col = 5, main = "Class-Metabolomics Correlations", filename = "plots/class_metabolomics_corr.png")

## Order Correlation matrices

# Order-Overall Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_order), as.matrix(df_sample_data_overall), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 7, main = "Order-Overall Correlations", filename = "plots/order_overall_corr.png")

# Order-Food Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_order), as.matrix(df_sample_data_food), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 7, main = "Order-Food Correlations", filename = "plots/order_food_corr.png")

# Order-Nutrients Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_order), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 7, fontsize_col = 4, main = "Order-Nutrients Correlations", filename = "plots/order_nutrients_corr.png")

# Order-Metabolomics Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_order), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 7, fontsize_col = 5, main = "Order-Metabolomics Correlations", filename = "plots/order_metabolomics_corr.png")

## Family Correlation matrices

# Family-Overall Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_family), as.matrix(df_sample_data_overall), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 3, main = "Family-Overall Correlations", filename = "plots/family_overall_corr.png")

# Family-Food Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_family), as.matrix(df_sample_data_food), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 3, main = "Family-Food Correlations", filename = "plots/family_food_corr.png")

# Family-Nutrients Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_family), as.matrix(df_sample_data_nutrients), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 3, fontsize_col = 4, main = "Family-Nutrients Correlations", filename = "plots/family_nutrients_corr.png")

# Family-Metabolomics Correlation matrix heatmap
cor_matrix <- cor(as.matrix(df_otu_table_family), as.matrix(df_sample_data_metabolomics), method = "spearman")
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, breaks = seq(-0.25, 0.25, length.out = length(default_palette) + 1), fontsize_row = 3, fontsize_col = 5, main = "Family-Metabolomics Correlations", filename = "plots/family_metabolomics_corr.png")

### PCA Analysis

## Phylum PCA Analysis

otu_table <- df_otu_table_phylum %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Phylum", values_to = "Abundance")
otu_table$Phylum[otu_table$Abundance < 0.10] <- "Other"
sample_data <- cbind(df_sample_data_overall, df_sample_data_food, df_sample_data_nutrients, df_sample_data_metabolomics)
num_sample_data <- sample_data[, sapply(sample_data, is.numeric)]  # Keep numeric columns
pca_result <- prcomp(num_sample_data, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(sample_data)
pca_with_otu <- left_join(otu_table, pca_df, by = "Sample")

explained_variance <- summary(pca_result)$importance[2, 1:2] * 100
x_label <- paste0("PC1 (", round(explained_variance[1], 2), "% Variance)")
y_label <- paste0("PC2 (", round(explained_variance[2], 2), "% Variance)")

ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2), color = "gray", size = 0.5) +  # Base points
  geom_arc_bar(
    data = pca_with_otu,
    aes(
      x0 = PC1, y0 = PC2, r0 = 0, r = 0.3,  # Radius controls pie chart size
      amount = Abundance, fill = Phylum
    ),
    stat = "pie"
  ) + coord_fixed() +
  labs(
    title = "PCA with Phylum Abundance Pie Charts",
    x = x_label,
    y = y_label,
    fill = "Phylum"
  ) +
  theme_minimal()
ggsave("plots/phylum_PCA.png")

## Class PCA Analysis

otu_table <- df_otu_table_class %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Class", values_to = "Abundance")
otu_table$Class[otu_table$Abundance < 0.10] <- "Other"
sample_data <- cbind(df_sample_data_overall, df_sample_data_food, df_sample_data_nutrients, df_sample_data_metabolomics)
num_sample_data <- sample_data[, sapply(sample_data, is.numeric)]  # Keep numeric columns
pca_result <- prcomp(num_sample_data, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(sample_data)
pca_with_otu <- left_join(otu_table, pca_df, by = "Sample")

explained_variance <- summary(pca_result)$importance[2, 1:2] * 100
x_label <- paste0("PC1 (", round(explained_variance[1], 2), "% Variance)")
y_label <- paste0("PC2 (", round(explained_variance[2], 2), "% Variance)")

ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2), color = "gray", size = 0.5) +  # Base points
  geom_arc_bar(
    data = pca_with_otu,
    aes(
      x0 = PC1, y0 = PC2, r0 = 0, r = 0.3,  # Radius controls pie chart size
      amount = Abundance, fill = Class
    ),
    stat = "pie"
  ) + coord_fixed() +
  labs(
    title = "PCA with Class Abundance Pie Charts",
    x = x_label,
    y = y_label,
    fill = "Class"
  ) +
  theme_minimal()
ggsave("plots/class_PCA.png")

## Order PCA Analysis

otu_table <- df_otu_table_order %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Order", values_to = "Abundance")
otu_table$Order[otu_table$Abundance < 0.10] <- "Other"
sample_data <- cbind(df_sample_data_overall, df_sample_data_food, df_sample_data_nutrients, df_sample_data_metabolomics)
num_sample_data <- sample_data[, sapply(sample_data, is.numeric)]  # Keep numeric columns
pca_result <- prcomp(num_sample_data, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(sample_data)
pca_with_otu <- left_join(otu_table, pca_df, by = "Sample")

explained_variance <- summary(pca_result)$importance[2, 1:2] * 100
x_label <- paste0("PC1 (", round(explained_variance[1], 2), "% Variance)")
y_label <- paste0("PC2 (", round(explained_variance[2], 2), "% Variance)")

ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2), color = "gray", size = 0.5) +  # Base points
  geom_arc_bar(
    data = pca_with_otu,
    aes(
      x0 = PC1, y0 = PC2, r0 = 0, r = 0.3,  # Radius controls pie chart size
      amount = Abundance, fill = Order
    ),
    stat = "pie"
  ) + coord_fixed() +
  labs(
    title = "PCA with Order Abundance Pie Charts",
    x = x_label,
    y = y_label,
    fill = "Order"
  ) +
  theme_minimal()
ggsave("plots/order_PCA.png")

## Family PCA Analysis

otu_table <- df_otu_table_family %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Family", values_to = "Abundance")
otu_table$Family[otu_table$Abundance < 0.10] <- "Other"
sample_data <- cbind(df_sample_data_overall, df_sample_data_food, df_sample_data_nutrients, df_sample_data_metabolomics)
num_sample_data <- sample_data[, sapply(sample_data, is.numeric)]  # Keep numeric columns
pca_result <- prcomp(num_sample_data, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(sample_data)
pca_with_otu <- left_join(otu_table, pca_df, by = "Sample")

explained_variance <- summary(pca_result)$importance[2, 1:2] * 100
x_label <- paste0("PC1 (", round(explained_variance[1], 2), "% Variance)")
y_label <- paste0("PC2 (", round(explained_variance[2], 2), "% Variance)")

ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2), color = "gray", size = 0.5) +  # Base points
  geom_arc_bar(
    data = pca_with_otu,
    aes(
      x0 = PC1, y0 = PC2, r0 = 0, r = 0.3,  # Radius controls pie chart size
      amount = Abundance, fill = Family
    ),
    stat = "pie"
  ) + coord_fixed() +
  labs(
    title = "PCA with Family Abundance Pie Charts",
    x = x_label,
    y = y_label,
    fill = "Family"
  ) +
  theme_minimal()
ggsave("plots/family_PCA.png")

# Features importance

loadings <- as.data.frame(pca_result$rotation)
loadings_squared <- loadings^2
loadings_squared$Feature <- rownames(loadings)
loadings_long <- tidyr::pivot_longer(loadings_squared, cols = starts_with("PC"), names_to = "PC", values_to = "Importance")

# Get top 50 features for PC1
top_pc1 <- loadings_long %>%
  filter(PC == "PC1") %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 50)

# Get top 50 features for PC2
top_pc2 <- loadings_long %>%
  filter(PC == "PC2") %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 50)

# Plot for PC1
ggplot(top_pc1, aes(x = reorder(Feature, -Importance), y = Importance, fill = Feature)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(
    title = "Top 50 Features for PC1",
    x = "Feature",
    y = "Importance (Squared Loadings)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("plots/importances_PC1.png")

# Plot for PC2
ggplot(top_pc2, aes(x = reorder(Feature, -Importance), y = Importance, fill = Feature)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(
    title = "Top 50 Features for PC2",
    x = "Feature",
    y = "Importance (Squared Loadings)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("plots/importances_PC2.png")