#############################################################################
#Title: "Rhizosphere Assembly - Heat map of Host-Selected ASVs"
#Author: "Ari Fina Bintarti"
#Date: "27-09-2026"
#############################################################################
# Selected ASVs = ASVs that are selected by Indicator Species analysis (IndVal.g) as well as Above Prediction by the Sloan Neutral Model analysis by stage

library(dplyr)
library(circlize)
library(grid)
library(ComplexHeatmap)

# Read the updated Uniform Taxonomy SILVA 144
new.tax.unif <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/uniform_taxonomy_edit.csv", row.names = 1)
new.tax.unif <- rownames_to_column(new.tax.unif, var = "asvID")
# Read the updated Weighted Taxonomy SILVA 144 (Weighted = based on specific habitat/sampling location, here I used the weighted taxonomy of Plant-Rhizosphere Soil = https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/weighted/plant-rhizosphere)
new.tax.weight <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/weight_taxonomy_edit.csv", row.names = 1)
new.tax.weight <- rownames_to_column(new.tax.weight, var = "asvID")

# I. RHIZOPLANE

# Read the rhizoplane indicator stage ASVs data which also contain all information from the other non-indicator stages
rp.indic.above.all.ed <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizoplane_indic_above_all.csv", row.names = 1)
# Edit the taxonomy columns
rp.indic.above.all.ed <- rp.indic.above.all.ed %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__")) %>%
  mutate(Weight.Class = str_remove(Weight.Class, "^c__")) %>%
  mutate(Weight.Order = str_remove(Weight.Order, "^o__")) %>%
  mutate(Weight.Family = str_remove(Weight.Family, "^f__")) %>%
  mutate(Weight.Genus = str_remove(Weight.Genus, "^g__")) %>%
  mutate(Weight.Species = str_remove(Weight.Species, "^s__")) 
# Make a new column with combination between Weight.Family and ASV ID for plotting
rp.indic.above.all.edit <- rp.indic.above.all.ed %>%
  mutate(Weight.Family_asvID = paste0(case_when(
    !is.na(Weight.Family) & Weight.Family != "" & Weight.Family != "--" & Weight.Family != "Incertae_Sedis" ~ Weight.Family,
    !is.na(Weight.Order) & Weight.Order != "" & Weight.Order != "--" & Weight.Order != "Incertae_Sedis" ~ Weight.Order,
    !is.na(Weight.Class) & Weight.Class != "" & Weight.Class != "--" & Weight.Class != "Incertae_Sedis" ~ Weight.Class,
    !is.na(Weight.Phylum) & Weight.Phylum != "" & Weight.Phylum != "--" & Weight.Phylum != "Incertae_Sedis" ~ Weight.Phylum, 
    TRUE ~ NA_character_), "-", substr(asvID, 1, 4)))
# Make some edits
rp.indic.above.all.edit[] <- lapply(rp.indic.above.all.edit, function(x) {
  if (is.character(x)) factor(x) else x
}) 
rp.indic.above.all.edit$Stage <- factor(rp.indic.above.all.edit$Stage,
                                      levels = c("Seedling","Early Vegetative",
                                                 "Late Vegetative","Early Reproductive",
                                                 "Late Reproductive"))
# Define stage order
stage_levels <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
# Clean important variables
rp.indic.above.all.edit <- rp.indic.above.all.edit %>%
  dplyr::mutate(Weight.Family_asvID = as.character(Weight.Family_asvID),
         Stage = as.character(Stage),
         Weight.Phylum = as.character(Weight.Phylum),
         stage_indicator = as.logical(stage_indicator))
# Check indicator stages
indicator_check <- rp.indic.above.all.edit %>%
  filter(stage_indicator == TRUE) %>%
  select(Weight.Family_asvID, Stage) %>%
  distinct() %>%
  group_by(Weight.Family_asvID) %>%
  summarise(n_indicator = n(), 
            indicator_stage = paste(Stage, collapse = "; "),
            .groups = "drop")
print(indicator_check)
# Check ASVs that do NOT have exactly one indicator stage
print(indicator_check %>%
    filter(n_indicator != 1))
# Create indicator-stage table
rp.indicator_stage <- rp.indic.above.all.edit %>%
  filter(stage_indicator == TRUE) %>%
  select(Weight.Family_asvID, Stage) %>%
  distinct() %>%
  dplyr::rename(indicator_stage = Stage)
# Convert to factor with desired order
rp.indicator_stage$indicator_stage <- factor(
  rp.indicator_stage$indicator_stage,
  levels = stage_levels)
# Prepare heat map data
rp.heatmap_data <- rp.indic.above.all.edit %>%
  select(Weight.Family_asvID, Stage, mean_relabund_percent, Weight.Phylum) %>%
  distinct() %>%
  left_join(rp.indicator_stage, by = "Weight.Family_asvID")
# Check ASVs without indicator stage
missing_indicator <- rp.heatmap_data %>%
  distinct(Weight.Family_asvID, indicator_stage) %>%
  filter(is.na(indicator_stage))
print(missing_indicator)
# Format stage for heat map
rp.heatmap_data <- rp.heatmap_data %>%
  mutate(Stage = factor(Stage,
                        levels = stage_levels,
                        labels = c("Seedling", "Early Vegetative", "Late Vegetative", 
                                   "Early Reproductive", "Late Reproductive")))
# Create ASV x Mean Relative Abundance matrix
rp.heatmap_matrix <- rp.heatmap_data %>%
  select(Weight.Family_asvID, Stage, mean_relabund_percent) %>%
  distinct() %>%
  pivot_wider(names_from = Stage,
              values_from = mean_relabund_percent, values_fill = 0) %>%
  as.data.frame()
# ASV IDs as row names
rownames(rp.heatmap_matrix) <- rp.heatmap_matrix$Weight.Family_asvID
rp.heatmap_matrix$Weight.Family_asvID <- NULL
# Convert to matrix
rp.heatmap_matrix <- as.matrix(rp.heatmap_matrix)
# Ensure correct stage column order
stage_order <- c("Seedling", "Early Vegetative", "Late Vegetative",
                 "Early Reproductive", "Late Reproductive")
rp.heatmap_matrix <- rp.heatmap_matrix[, stage_order, drop = FALSE]
# Remove ASV with zero variance
# Z-score cannot be calculated if an ASV has exactly the same abundance in all five stages.
asv_sd <- apply(rp.heatmap_matrix, 1, sd, na.rm = TRUE)
keep_asv <- !is.na(asv_sd) &
  asv_sd > 0
rp.heatmap_matrix <- rp.heatmap_matrix[
  keep_asv,,drop = FALSE]
# Calculate row-wise Z score
# Z-score is calculated across the 5 stages for EACH ASV.
rp.heatmap_matrix_z <- t(scale(t(rp.heatmap_matrix)))
# Restore row and column names
rownames(rp.heatmap_matrix_z) <- rownames(rp.heatmap_matrix)
colnames(rp.heatmap_matrix_z) <- colnames(rp.heatmap_matrix)
# Create row information
rp.row_info <- rp.heatmap_data %>%
  select(Weight.Family_asvID, indicator_stage, Weight.Phylum) %>%
  distinct() %>%
  mutate(Weight.Family_asvID = as.character(Weight.Family_asvID),
    indicator_stage = as.character(indicator_stage),
    Weight.Phylum = as.character(Weight.Phylum)) %>%
  filter(Weight.Family_asvID %in%
      rownames(rp.heatmap_matrix_z))
# Check for duplicate ASVs
duplicate_asv <- rp.row_info %>%
  count(Weight.Family_asvID) %>%
  filter(n > 1)
print(duplicate_asv)
# Define indicator stage order
indicator_order <- c("Seedling","Early Vegetative","Late Vegetative",
                     "Early Reproductive","Late Reproductive")
# Order ASVs by:
# 1. Indicator stage
# 2. Phylum according to phylum_cols order
# 3. ASV ID

# Get phylum order from your legend
phylum_cols <- readRDS("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/phylum_cols.rds")
# Get the Phylum order 
phylum_order <- names(phylum_cols)
# Make Phylum a factor using the SAME order as phylum_cols
rp.row_info <- rp.row_info %>%
  mutate(indicator_stage = factor(indicator_stage,
                                  levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum,
                                     levels = phylum_order)) %>%
  arrange(indicator_stage, Weight.Phylum, Weight.Family_asvID)
# Create the definitive ASV order
asv_order <- as.character(rp.row_info$Weight.Family_asvID)
#rp.row_labels <- asv_order %>%
  #sub("^(.*)-([^-]+)$", "\\2-\\1", .)
# Apply exact same ASV order to Z-score matrix
rp.heatmap_matrix_z <- rp.heatmap_matrix_z[
  asv_order,, drop = FALSE]
# Check synchronization
check_sync <- identical(rownames(rp.heatmap_matrix_z), 
                        as.character(rp.row_info$Weight.Family_asvID))
print(check_sync)
# Add a check table
rp.check_order <- data.frame(
  ASV = rownames(rp.heatmap_matrix_z),
  Indicator_stage = as.character(rp.row_info$indicator_stage),
  Phylum = as.character(rp.row_info$Weight.Phylum),
  stringsAsFactors = FALSE)
print(head(rp.check_order, 20))
# Phylum colors and sort by alphabetical order
phylum_used <- sort(unique(as.character(rp.row_info$Weight.Phylum)))
# Match the colors to the actual Phyla used
phylum_cols_used <- phylum_cols[phylum_used]
# Restore names
names(phylum_cols_used) <- phylum_used
# Check for missing colors
if (any(is.na(phylum_cols_used))) {
   missing_phyla <- names(phylum_cols_used)[is.na(phylum_cols_used)]
  stop(paste("Missing Phylum colors:",
             paste(missing_phyla, collapse = ", ")))
}
# Indicator stage colors
indicator_stage_cols <- c(
  "Seedling" = "lightgrey",
  "Early Vegetative" = "lightgrey",
  "Late Vegetative" = "lightgrey",
  "Early Reproductive" = "lightgrey",
  "Late Reproductive" = "lightgrey")
hcl.colors(5, palette = "ag_GrnYl")
# Z-Score color scale
rp.z_col <- colorRamp2(c(-2, 0, 2),
  c("#2166AC","white","#B2182B"))
# Create row-annotation
rp.ha <- rowAnnotation(Phylum = rp.row_info$Weight.Phylum,
                       `Indicator stage` = rp.row_info$indicator_stage,
                       col = list(Phylum = phylum_cols_used, `Indicator stage` = indicator_stage_cols),
                       show_annotation_name = c(TRUE, FALSE),
                       annotation_name_gp = gpar(fontsize = 15), 
                       gap = unit(0.8, "mm"))
#rp.ha <- rowAnnotation(
  #Phylum = rp.row_info$Weight.Phylum, 
  #col = list(Phylum = phylum_cols_used),
  #show_annotation_name = TRUE,
  #annotation_name_gp = gpar(fontsize = 15))

# Create Heatmap
cn <- c("Seedling", "Early\nVegetative", "Late\nVegetative",
  "Early\nReproductive", "Late\nReproductive")
#ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
ht_opt$TITLE_PADDING= unit(6,"mm") # general safety cushion
ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                  # space between annotation and heatmap
rp.ht <- Heatmap(rp.heatmap_matrix_z, 
                 name = "Z-score",
                 col = rp.z_col,
                 # Phylum annotation → RIGHT
                 right_annotation = rp.ha,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 # Indicator stage labels → RIGHT
                 row_split = rp.row_info$indicator_stage,
                 row_title_side = "right",
                 # ASV names → LEFT
                 show_row_names = TRUE,
                 row_names_side = "left",
                 row_labels = asv_order,#rp.row_labels,
                 show_column_names = FALSE,
                 row_names_max_width = unit(7, "cm"),
                 row_names_gp = gpar(fontsize = 9.6),
                 row_title_gp = gpar(col = "black", fontsize = 16, fontface = "bold"), #fill = "grey95"
                 gap = unit(3, "mm"), 
                 border = TRUE,
                 border_gp = gpar(col = "#1B1B1B", lwd = 0.8),
                 width = unit(65, "mm"),
                 column_title = "A. Rhizoplane",
                 column_title_gp = gpar(fontsize = 28),
                 bottom_annotation = HeatmapAnnotation(
                   text = anno_text(cn, rot = 90, location = unit(0.9, "npc"), just = "right",
                                    gp = gpar(fontsize = 14)), 
                   annotation_height = unit(30, "mm")),
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 14),
                                             labels_gp = gpar(fontsize = 12)))
# Create Z-score legend
#lgd_z <- Legend(title = "Z-score", 
                #col_fun = rp.z_col, 
                #title_gp = gpar(fontsize = 14),
                #labels_gp = gpar(fontsize = 12))
# Create indicator stage legend
#lgd_stage <- Legend(title = "Indicator stage",
                   #at = indicator_order,
                   #legend_gp = gpar(fill = indicator_stage_cols[indicator_order]),
                   #title_gp = gpar(fontsize = 14),
                   #labels_gp = gpar(fontsize = 12))
# Create Phylum legend
#lgd_phylum <- Legend(title = "Phylum",
                     #at = names(phylum_cols_used),
                     #legend_gp = gpar(fill = phylum_cols_used),
                     #title_gp = gpar(fontsize = 14),
                     #labels_gp = gpar(fontsize = 12),
                     #grid_height = unit(4, "mm"),
                     #grid_width = unit(4, "mm"))
# Pack all 3 legends vertically
#all_legends <- packLegend(lgd_z, lgd_phylum, 
                          #direction = "vertical",
                          #gap = unit(8, "mm"))
# Save the plot
#setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/")
#tiff("Rhizoplane_selectASV_Zscore_heatmap.tiff",
     #width = 13, 
     #height = 20, 
     #units = "in", 
     #res = 300, 
     #compression = "lzw")
#draw(rp.ht,
     #show_heatmap_legend = FALSE,
     #show_annotation_legend = FALSE)
#dev.off()


# II. RHIZOSPHERE

# Read the edited data
rs.indic.above.all.ed <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizosphere_indic_above_all.csv", row.names = 1)
# Edit the taxonomy columns
rs.indic.above.all.ed <- rs.indic.above.all.ed %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__")) %>%
  mutate(Weight.Class = str_remove(Weight.Class, "^c__")) %>%
  mutate(Weight.Order = str_remove(Weight.Order, "^o__")) %>%
  mutate(Weight.Family = str_remove(Weight.Family, "^f__")) %>%
  mutate(Weight.Genus = str_remove(Weight.Genus, "^g__")) %>%
  mutate(Weight.Species = str_remove(Weight.Species, "^s__")) 
# Make a new column with combination between Weight.Family and ASV ID for plotting
rs.indic.above.all.edit <- rs.indic.above.all.ed %>%
  mutate(Weight.Family_asvID = paste0(case_when(
    !is.na(Weight.Family) & Weight.Family != "" & Weight.Family != "--" & Weight.Family != "Incertae_Sedis" ~ Weight.Family,
    !is.na(Weight.Order) & Weight.Order != "" & Weight.Order != "--" & Weight.Order != "Incertae_Sedis" ~ Weight.Order,
    !is.na(Weight.Class) & Weight.Class != "" & Weight.Class != "--" & Weight.Class != "Incertae_Sedis" ~ Weight.Class,
    !is.na(Weight.Phylum) & Weight.Phylum != "" & Weight.Phylum != "--" & Weight.Phylum != "Incertae_Sedis" ~ Weight.Phylum, 
    TRUE ~ NA_character_), "-", substr(asvID, 1, 4)))
# Make some edits
rs.indic.above.all.edit[] <- lapply(rs.indic.above.all.edit, function(x) {
  if (is.character(x)) factor(x) else x
}) 
rs.indic.above.all.edit$Stage <- factor(rs.indic.above.all.edit$Stage,
                                      levels = c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive"))
# Define stage order
stage_levels <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
# Clean important variables
rs.indic.above.all.edit <- rs.indic.above.all.edit %>%
  mutate(Weight.Family_asvID = as.character(Weight.Family_asvID),
         Stage = as.character(Stage),
         Weight.Phylum = as.character(Weight.Phylum),
         stage_indicator = as.logical(stage_indicator))
# Check indicator stages
indicator_check <- rs.indic.above.all.edit %>%
  filter(stage_indicator == TRUE) %>%
  select(Weight.Family_asvID, Stage) %>%
  distinct() %>%
  group_by(Weight.Family_asvID) %>%
  summarise(n_indicator = n(), 
            indicator_stage = paste(Stage, collapse = "; "),
            .groups = "drop")
print(indicator_check)
# Check ASVs that do NOT have exactly one indicator stage
print(indicator_check %>%
        filter(n_indicator != 1))
# Create indicator-stage table
rs.indicator_stage <- rs.indic.above.all.edit %>%
  filter(stage_indicator == TRUE) %>%
  select(Weight.Family_asvID, Stage) %>%
  distinct() %>%
  dplyr::rename(indicator_stage = Stage)
# Convert to factor with desired order
rs.indicator_stage$indicator_stage <- factor(
  rs.indicator_stage$indicator_stage,
  levels = stage_levels)
# Prepare heat map data
rs.heatmap_data <- rs.indic.above.all.edit %>%
  select(Weight.Family_asvID, Stage, mean_relabund_percent, Weight.Phylum) %>%
  distinct() %>%
  left_join(rs.indicator_stage, by = "Weight.Family_asvID")
# Format stage for heat map
rs.heatmap_data <- rs.heatmap_data %>%
  mutate(Stage = factor(Stage,
                        levels = stage_levels,
                        labels = c("Seedling", "Early Vegetative", "Late Vegetative", 
                                   "Early Reproductive", "Late Reproductive")))
# Create ASV x Mean Relative Abundance matrix
rs.heatmap_matrix <- rs.heatmap_data %>%
  select(Weight.Family_asvID, Stage, mean_relabund_percent) %>%
  distinct() %>%
  pivot_wider(names_from = Stage,
              values_from = mean_relabund_percent, values_fill = 0) %>%
  as.data.frame()
# ASV IDs as row names
rownames(rs.heatmap_matrix) <- rs.heatmap_matrix$Weight.Family_asvID
rs.heatmap_matrix$Weight.Family_asvID <- NULL
# Convert to matrix
rs.heatmap_matrix <- as.matrix(rs.heatmap_matrix)
# Ensure correct stage column order
stage_order <- c("Seedling", "Early Vegetative", "Late Vegetative",
                 "Early Reproductive", "Late Reproductive")
rs.heatmap_matrix <- rs.heatmap_matrix[, stage_order, drop = FALSE]
# Remove ASV with zero variance
# Z-score cannot be calculated if an ASV has exactly the same abundance in all five stages.
asv_sd <- apply(rs.heatmap_matrix, 1, sd, na.rm = TRUE)
keep_asv <- !is.na(asv_sd) &
  asv_sd > 0
rs.heatmap_matrix <- rs.heatmap_matrix[
  keep_asv,,drop = FALSE]
# Calculate row-wise Z score
# Z-score is calculated across the 5 stages for EACH ASV.
rs.heatmap_matrix_z <- t(scale(t(rs.heatmap_matrix)))
# Restore row and column names
rownames(rs.heatmap_matrix_z) <- rownames(rs.heatmap_matrix)
colnames(rs.heatmap_matrix_z) <- colnames(rs.heatmap_matrix)
# Create row information
rs.row_info <- rs.heatmap_data %>%
  select(Weight.Family_asvID, indicator_stage, Weight.Phylum) %>%
  distinct() %>%
  mutate(Weight.Family_asvID = as.character(Weight.Family_asvID),
         indicator_stage = as.character(indicator_stage),
         Weight.Phylum = as.character(Weight.Phylum)) %>%
  filter(Weight.Family_asvID %in%
           rownames(rs.heatmap_matrix_z))
# Define indicator stage order
indicator_order <- c("Seedling","Early Vegetative","Late Vegetative",
                     "Early Reproductive","Late Reproductive")
# Order ASVs by:
# 1. Indicator stage
# 2. Phylum according to phylum_cols order
# 3. ASV ID

# Get the Phylum order 
phylum_order <- names(phylum_cols)
# Make Phylum a factor using the SAME order as phylum_cols
rs.row_info <- rs.row_info %>%
  mutate(indicator_stage = factor(indicator_stage,
                                  levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum,
                                     levels = phylum_order)) %>%
  arrange(indicator_stage, Weight.Phylum, Weight.Family_asvID)
# Create the definitive ASV order
asv_order <- as.character(rs.row_info$Weight.Family_asvID)
# Labels displayed in the heatmap only
#rs.row_labels <- asv_order %>%
  #sub("^(.*)-([^-]+)$", "\\2-\\1", .)
# Apply exact same ASV order to Z-score matrix
rs.heatmap_matrix_z <- rs.heatmap_matrix_z[
  asv_order,, drop = FALSE]
# Phylum colors and sort by alphabetical order
phylum_used <- sort(unique(as.character(rs.row_info$Weight.Phylum)))
# Match the colors to the actual Phyla used
phylum_cols_used <- phylum_cols[phylum_used]
# Restore names
names(phylum_cols_used) <- phylum_used
# Indicator stage colors
indicator_stage_cols <- c(
  "Seedling" = "lightgrey",
  "Early Vegetative" = "lightgrey",
  "Late Vegetative" = "lightgrey",
  "Early Reproductive" = "lightgrey",
  "Late Reproductive" = "lightgrey")
hcl.colors(5, palette = "ag_GrnYl")
# Z-Score color scale
rs.z_col <- colorRamp2(c(-2, 0, 2),
                       c("#2166AC","white","#B2182B"))
# Create row-annotation
rs.ha <- rowAnnotation(Phylum = rs.row_info$Weight.Phylum,
                       `Indicator stage` = rs.row_info$indicator_stage,
                        col = list(Phylum = phylum_cols_used,
                                  `Indicator stage` = indicator_stage_cols),
                       show_annotation_name = c(TRUE, FALSE),
                       annotation_name_gp = gpar(fontsize = 15), gap = unit(0.5, "mm"))
#rs.ha <- rowAnnotation(
  #Phylum = rs.row_info$Weight.Phylum, 
  #col = list(Phylum = phylum_cols_used),
  #show_annotation_name = TRUE,
  #annotation_name_gp = gpar(fontsize = 15))

# Create Heatmap
cn <- c("Seedling", "Early\nVegetative", "Late\nVegetative",
        "Early\nReproductive", "Late\nReproductive")
rs.ht <- Heatmap(rs.heatmap_matrix_z, 
                 name = "Z-score",
                 col = rs.z_col,
                 # Phylum annotation → RIGHT
                 right_annotation = rs.ha,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 # Indicator stage labels → RIGHT
                 row_split = rs.row_info$indicator_stage,
                 row_title_side = "right",
                 # ASV names → LEFT
                 show_row_names = TRUE,
                 row_names_side = "left",
                 row_labels = asv_order,#rs.row_labels,
                 show_column_names = FALSE,
                 row_names_max_width = unit(7, "cm"),
                 row_names_gp = gpar(fontsize = 9.6),
                 row_title_gp = gpar(col =  "black", fontsize = 16, fontface = "bold"), #fill = "grey95",
                 gap = unit(3, "mm"), 
                 border = TRUE,
                 border_gp = gpar(col = "#1B1B1B", lwd = 0.8),
                 width = unit(65, "mm"),
                 column_title = "B. Rhizosphere",
                 column_title_gp = gpar(fontsize = 28),
                 bottom_annotation = HeatmapAnnotation(
                   text = anno_text(cn, 
                                    rot = 90, 
                                    location = unit(0.9, "npc"), 
                                    just = "right",
                                    gp = gpar(fontsize = 14)),
                   annotation_height = unit(30, "mm")),
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 14),
                                             labels_gp = gpar(fontsize = 12)))
# Create Z-score legend
#lgd_z <- Legend(title = "Z-score", 
                #col_fun = rs.z_col, 
                #title_gp = gpar(fontsize = 14),
                #labels_gp = gpar(fontsize = 12))
# Create indicator stage legend
#lgd_stage <- Legend(title = "Indicator stage",
                    #at = indicator_order,
                    #legend_gp = gpar(fill = indicator_stage_cols[indicator_order]),
                    #title_gp = gpar(fontsize = 14),
                    #labels_gp = gpar(fontsize = 12))
# Create Phylum legend
#lgd_phylum <- Legend(title = "Phylum",
                     #at = names(phylum_cols_used),
                     #legend_gp = gpar(fill = phylum_cols_used),
                     #title_gp = gpar(fontsize = 14),
                     #labels_gp = gpar(fontsize = 12),
                     #grid_height = unit(4, "mm"),
                     #grid_width = unit(4, "mm"))
# Pack all 3 legends vertically
#all_legends <- packLegend(lgd_z, lgd_phylum, 
                          #direction = "vertical",
                          #gap = unit(8, "mm"))
# Save the plot
#setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/")
#tiff("Rhizosphere_selectASV_Zscore_heatmap.tiff",
     #width = 13, 
     #height = 20, 
     #units = "in", 
     #res = 300, 
     #compression = "lzw")
#draw(rs.ht,
     #show_heatmap_legend = FALSE,
     #show_annotation_legend = FALSE)
#dev.off()

#=========================================================#
# Combine Both Figures
#=========================================================#

# Create Common Phylum Legend
lgd_phylum_common <- Legend(title = "Phylum", at = names(phylum_cols),
  legend_gp = gpar(fill = phylum_cols),
  title_gp = gpar(fontsize = 15),
  labels_gp = gpar(fontsize = 14),
  row_gap = unit(2, "mm"),
  title_gap = unit(4, "mm"))
# Create Common Z-score legend
lgd_z_common <- Legend(title = "Z-score",col_fun = rp.z_col, 
                       title_gp = gpar(fontsize = 15),
                       labels_gp = gpar(fontsize = 14),
                       title_gap = unit(4, "mm"))
# Combine 2 legends vertically
all_legends_common <- packLegend(lgd_z_common, lgd_phylum_common,
                                 direction = "vertical", gap = unit(8, "mm"))
# Grab the legend
hm.legend_grob <- grid.grabExpr(draw(all_legends_common))
# Grab each heatmap separately
rp1 <- grid.grabExpr(draw(rp.ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE))
rs2 <- grid.grabExpr(draw(rs.ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE))
# Combine heatmap and common legend
#hm_combined <- plot_grid(rp1, rs2, hm.legend_grob, ncol = 3, rel_widths = c(1, 0.9, 0.18), align = "h")
hm_combined <- ggdraw() +
  draw_plot(rp1, x = -0.06,    y = 0, width = 0.49, height = 1) +
  draw_plot(rs2, x = 0.36, y = 0, width = 0.49, height = 1) +
  draw_plot(hm.legend_grob, x = 0.83, y = 0.05, width = 0.14, height = 0.9)
# Calculate height for saving plot
rp_height <- length(unique(rp.indic.above.all.edit$Weight.Family_asvID)) * 0.125
rs_height <- length(unique(rs.indic.above.all.edit$Weight.Family_asvID)) * 0.125
combined_height <- max(rp_height, rs_height)
# Save Figure
setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/Tidy_R_Figures/")
#tiff("FigXX_ZScore_HeatMap_byStage.tiff",
  #width = 16,
  #height = combined_height, 
  #units = "in",
  #res = 300,
  #type = "cairo", 
  #compression = "lzw")
#print(hm_combined)
#dev.off()



