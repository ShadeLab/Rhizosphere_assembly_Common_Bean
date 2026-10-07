#############################################################################
#Title: "Rhizosphere Assembly - Bubble Plot of Host-Selected ASVs"
#Author: "Ari Fina Bintarti"
#Date: "26-09-2026"
#############################################################################
# Host-Selected ASVs = Overlapped/shared ASVs between Indicator Species analysis (IndVal.g) and Sloan Neutral Model = Above Prediction

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(scales)

# Read the updated Uniform Taxonomy SILVA 144
new.tax.unif <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/Taxonomy_SILVA_144/uniform_taxonomy_edit.csv", row.names = 1)
new.tax.unif <- rownames_to_column(new.tax.unif, var = "asvID")
# Read the updated Weighted Taxonomy SILVA 144 (Weighted = based on specific habitat/sampling location, here I used the weighted taxonomy of Plant-Rhizosphere Soil = https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/weighted/plant-rhizosphere)
new.tax.weight <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/Taxonomy_SILVA_144/weight_taxonomy_edit.csv", row.names = 1)
new.tax.weight <- rownames_to_column(new.tax.weight, var = "asvID")

#######################################################################
### Prepare Data Set for Plotting Bubble Plot of Host-Selected Taxa ###
#######################################################################

# I. Rhizoplane

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
# change stage level
rp.indic.above.all.edit$Stage <- factor(
  rp.indic.above.all.edit$Stage, levels = c("Seedling","Early Vegetative","Late Vegetative",
                                          "Early Reproductive","Late Reproductive"))
# Define stage order
stage_order <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
indicator_order <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
# Identify the indicator stage for each ASV
rp.indic.above.all.edit <- rp.indic.above.all.edit %>%
  group_by(Weight.Family_asvID) %>%
  mutate(indicator_stage = Stage[stage_indicator][1]) %>%
  ungroup()
# Make data for bubble plot
rp.bubble_data <- rp.indic.above.all.edit %>%
  mutate(Stage = factor(Stage, levels = stage_order),
         indicator_stage = factor(indicator_stage, levels = indicator_order))
# Get phylum order from your legend
phylum_cols <- readRDS("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/phylum_cols.rds")
phylum_order <- names(phylum_cols)
# Create ASV information table
rp.bubble.asv_info <- rp.bubble_data %>% 
  distinct(Weight.Family_asvID, indicator_stage, Weight.Phylum) %>%
  mutate(indicator_stage = factor(indicator_stage, levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum, levels = phylum_order))
# Order ASVs WITHIN each indicator-stage facet
rp.bubble.asv_order <- rp.bubble.asv_info %>% 
  arrange(indicator_stage, Weight.Phylum, Weight.Family_asvID) %>%
  pull(Weight.Family_asvID)
# Reverse because ggplot displays the first factor level at the bottom
rp.bubble.asv_order <- rev(rp.bubble.asv_order)
# Apply ASV order
rp.bubble_data <- rp.bubble_data %>%
  mutate(Weight.Family_asvID = factor(Weight.Family_asvID, levels = rp.bubble.asv_order),
         Stage = factor(Stage, levels = stage_order))
# Make the bubble Plot
rp.bubblePlot <- ggplot(
  rp.bubble_data, aes(x = Stage, y = Weight.Family_asvID)) +
  geom_point(data = ~ subset(.x, stage_indicator == TRUE),
             aes(size = mean_relabund_percent, fill = Weight.Phylum), shape = 21, colour = "black", alpha = 0.8) +
  geom_point(data = ~ subset(.x, stage_indicator == FALSE), aes(size = mean_relabund_percent),
             shape = 21, fill = "white",colour = "black",alpha = 0.6) +
  scale_size_continuous(name = "Mean Relative\nAbundance (%)", range = c(2, 8)) +
  scale_fill_manual(values = phylum_cols, breaks = names(phylum_cols), drop = TRUE) +
  facet_grid(indicator_stage ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y= NULL, fill="Phylum") +
  #scale_x_discrete(labels = c("Seedling" = "Seedling",
      #"Early Vegetative" = "Early\nVegetative",
      #"Late Vegetative" = "Late\nVegetative",
      #"Early Reproductive" = "Early\nReproductive",
      #"Late Reproductive" = "Late\nReproductive")) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14), # 14
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 13, face = "bold"),
        strip.background = element_rect(colour = "lightgrey"),
        legend.position = "none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        panel.border = element_rect(colour = "lightgrey", fill = NA,linewidth = 0.5),
        plot.margin = margin(10, 10, 1, 10)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 6)))
# Make the IndVal bar plot data
rp.indval_data <- rp.bubble_data %>%
  filter(stage_indicator == TRUE) %>%
  distinct(Weight.Family_asvID, indicator_stage, Weight.Phylum, IndVal.Stat) %>%
  filter(!is.na(IndVal.Stat)) %>%
  mutate(Weight.Family_asvID = factor(Weight.Family_asvID, levels = levels(rp.bubble_data$Weight.Family_asvID)),
         indicator_stage = factor(indicator_stage, levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum, levels = names(phylum_cols)))
# Check if the ASVs order of the bar plot data identical with the bubble data
identical(levels(rp.bubble_data$Weight.Family_asvID),
          levels(rp.indval_data$Weight.Family_asvID))
# Make the IndVal bar plot
rp.indvalPlot <- ggplot(rp.indval_data,
                        aes(x = IndVal.Stat, y = Weight.Family_asvID, fill = Weight.Phylum)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = phylum_cols, breaks = names(phylum_cols), drop = TRUE) +
  facet_grid(indicator_stage ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), labels = scales::label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "IndVal", y = NULL, fill = "Phylum") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15, margin = margin(t = -50)),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "lightgrey", fill = NA, linewidth = 0.5),
        panel.spacing.y = unit(2, "mm"),
        plot.margin = margin(10, -5, 1, 10))
# Combine both figures and Save
rp.combined <- rp.indvalPlot + rp.bubblePlot +
  plot_layout(widths = c(0.4, 1))
# Add the common title of Rhizoplane
rp.combined_w_title <- ggdraw() +
  draw_plot(plot_grid(ggdraw() + draw_label("A. Rhizoplane", x = 0, y = 0.98, size = 28, hjust=0, vjust=1),
                      rp.combined, ncol = 1, rel_heights = c(0.025, 1)), x = 0, y = 0.01, width = 0.99,height = 0.99)
# Save the plot
#setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/")
#ggsave("Rhizoplane_RA_IndValPlot.tiff",
       #rp.combined_w_title , device = "tiff",
       #width = 8, 
       #height = length(unique(rp.indic.above.all.edit$Weight.Family_asvID)) * 0.11,
       #units= "in", dpi = 300,
       #compression="lzw", bg= "white")


# II. Rhizosphere

# Read the rhizosphere indicator stage ASVs data which also contain all information from the other non-indicator stages
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
# change stage level
rs.indic.above.all.edit$Stage <- factor(
  rs.indic.above.all.edit$Stage, levels = c("Seedling","Early Vegetative","Late Vegetative",
                                            "Early Reproductive","Late Reproductive"))
# Define stage order
stage_order <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
indicator_order <- c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive")
# Identify the indicator stage for each ASV
rs.indic.above.all.edit <- rs.indic.above.all.edit %>%
  group_by(Weight.Family_asvID) %>%
  mutate(indicator_stage = Stage[stage_indicator][1]) %>%
  ungroup()
# Make data for bubble plot
rs.bubble_data <- rs.indic.above.all.edit %>%
  mutate(Stage = factor(Stage, levels = stage_order),
         indicator_stage = factor(indicator_stage, levels = indicator_order))
# Get phylum order from your legend
phylum_cols <- readRDS("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/phylum_cols.rds")
phylum_order <- names(phylum_cols)
# Create ASV information table
rs.bubble.asv_info <- rs.bubble_data %>% 
  distinct(Weight.Family_asvID, indicator_stage, Weight.Phylum) %>%
  mutate(indicator_stage = factor(indicator_stage, levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum, levels = phylum_order))
# Order ASVs WITHIN each indicator-stage facet
rs.bubble.asv_order <- rs.bubble.asv_info %>% 
  arrange(indicator_stage, Weight.Phylum, Weight.Family_asvID) %>%
  pull(Weight.Family_asvID)
# Reverse because ggplot displays the first factor level at the bottom
rs.bubble.asv_order <- rev(rs.bubble.asv_order)
# Apply ASV order
rs.bubble_data <- rs.bubble_data %>%
  mutate(Weight.Family_asvID = factor(Weight.Family_asvID, levels = rs.bubble.asv_order),
         Stage = factor(Stage, levels = stage_order))
# Make the bubble Plot
rs.bubblePlot <- ggplot(
  rs.bubble_data, aes(x = Stage, y = Weight.Family_asvID)) +
  geom_point(data = ~ subset(.x, stage_indicator == TRUE),
             aes(size = mean_relabund_percent, fill = Weight.Phylum), shape = 21, colour = "black", alpha = 0.8) +
  geom_point(data = ~ subset(.x, stage_indicator == FALSE), aes(size = mean_relabund_percent),
             shape = 21, fill = "white",colour = "black",alpha = 0.6) +
  scale_size_continuous(name = "Mean Relative\nAbundance (%)", range = c(2, 8)) +
  scale_fill_manual(values = phylum_cols, breaks = names(phylum_cols), drop = TRUE) +
  facet_grid(indicator_stage ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y= NULL, fill="Phylum") +
  #scale_x_discrete(labels = c("Seedling" = "Seedling",
                              #"Early Vegetative" = "Early\nVegetative",
                              #"Late Vegetative" = "Late\nVegetative",
                              #"Early Reproductive" = "Early\nReproductive",
                             # "Late Reproductive" = "Late\nReproductive")) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14), # 14
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(size = 13, face = "bold"),
        strip.background = element_rect(colour = "lightgrey"),
        legend.position = "none",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        panel.border = element_rect(colour = "lightgrey", fill = NA,linewidth = 0.5),
        plot.margin = margin(10, 10, 1, 10)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 6)))
# Make the IndVal bar plot data
rs.indval_data <- rs.bubble_data %>%
  filter(stage_indicator == TRUE) %>%
  distinct(Weight.Family_asvID, indicator_stage, Weight.Phylum, IndVal.Stat) %>%
  filter(!is.na(IndVal.Stat)) %>%
  mutate(Weight.Family_asvID = factor(Weight.Family_asvID, levels = levels(rs.bubble_data$Weight.Family_asvID)),
         indicator_stage = factor(indicator_stage, levels = indicator_order),
         Weight.Phylum = factor(Weight.Phylum, levels = names(phylum_cols)))
# Check if the ASVs order of the bar plot data identical with the bubble data
identical(levels(rs.bubble_data$Weight.Family_asvID),
          levels(rs.indval_data$Weight.Family_asvID))
# Make the IndVal bar plot
rs.indvalPlot <- ggplot(rs.indval_data,
                        aes(x = IndVal.Stat, y = Weight.Family_asvID, fill = Weight.Phylum)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = phylum_cols, breaks = names(phylum_cols), drop = TRUE) +
  facet_grid(indicator_stage ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), labels = scales::label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "IndVal", y = NULL, fill = "Phylum") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15, margin = margin(t = -50)),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "lightgrey", fill = NA, linewidth = 0.5),
        panel.spacing.y = unit(2, "mm"),
        plot.margin = margin(10, -5, 1, 10))
# Combine both figures and Save
rs.combined <- rs.indvalPlot + rs.bubblePlot +
  plot_layout(widths = c(0.4, 1))
# Add the common title of Rhizosphere
rs.combined_w_title <- ggdraw() +
  draw_plot(plot_grid(ggdraw() + draw_label("B. Rhizosphere", x = 0, y = 0.98, size = 28, hjust=0, vjust=1),
                      rs.combined, ncol = 1, rel_heights = c(0.025, 1)), x = 0, y = 0.01, width = 0.99,height = 0.99)
# Save the plot
#setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/09072026/")
#ggsave("Rhizosphere_RA_IndValPlot.tiff",
       #rs.combined_w_title , device = "tiff",
       #width = 8, 
       #height = length(unique(rs.indic.above.all.edit$Weight.Family_asvID)) * 0.11,
       #units= "in", dpi = 300,
       #compression="lzw", bg= "white")

# Common Phylum legend
phylum_legend_dummy <- ggplot(
  data.frame(Weight.Phylum = factor(names(phylum_cols), levels = names(phylum_cols)), 
             x = seq_along(phylum_cols), y = 1),
  aes(x = x, y = y, fill = Weight.Phylum)) +
  geom_point(shape = 21, size = 5.5, colour = "black") +
  scale_fill_manual(values = phylum_cols, breaks = names(phylum_cols), drop = FALSE) +
  labs(fill = "Phylum") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.justification = "left",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        legend.title.position = "top",
        legend.box.just = "left",
        legend.key.height = unit(0.3, "cm")) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, 
                             override.aes = list(shape = 21, size = 6)))

# Extract common legend Phylum
phylum_legend <- cowplot::get_legend(phylum_legend_dummy)

# Mean Relative Abundance Legend
# Combine both data
all_relabund <- c(rp.bubble_data$mean_relabund_percent, rs.bubble_data$mean_relabund_percent)
all_relabund <- all_relabund[is.finite(all_relabund)]
# Choose the break
size_breaks <- pretty(all_relabund, n = 5)
size_breaks <- size_breaks[size_breaks >= min(all_relabund) &
                             size_breaks <= max(all_relabund)]
size_breaks
# Make mean relative abundance dummy plot
size_legend_dummy <- ggplot(data.frame(x = seq_along(size_breaks), y = 1, 
                                       mean_relabund_percent = size_breaks),
                            aes(x = x, y = y, size = mean_relabund_percent)) +
  geom_point(shape = 21, fill = "white", colour = "black") +
  scale_size_continuous(name = "Mean Relative\nAbundance (%)", 
                        breaks = size_breaks, range = c(2, 8)) +
  theme_void() +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.box.just = "left")
# Extract common legend Mean Relative Abundance
size_legend <- cowplot::get_legend(size_legend_dummy)
# Calculate height for saving plot
rp_height <- length(unique(rp.indic.above.all.edit$Weight.Family_asvID)) * 0.12
rs_height <- length(unique(rs.indic.above.all.edit$Weight.Family_asvID)) * 0.12
# Combine the Rhizoplane and Rhizosphere bubble plots
rp.combined_w_title <- rp.combined_w_title +
  theme(plot.margin = margin(0, 10, 0, 0))
rs.combined_w_title <- rs.combined_w_title +
  theme(plot.margin = margin(0, 0, 0, 12))
combined.bubble.plot <- cowplot::plot_grid(rp.combined_w_title, rs.combined_w_title,
                                           ncol = 2, align = "none")
# Combine the two legends
legends_combined <- cowplot::plot_grid(size_legend, phylum_legend,  
                                       ncol = 1, rel_heights = c(0.45, 1))
# Final figure
combined.bubble.final <- cowplot::plot_grid(combined.bubble.plot, 
                                            cowplot::plot_grid(NULL, legends_combined, NULL,ncol = 1,
                                                               rel_heights = c(0.15, 0.70, 0.15)), 
                                            ncol = 2, rel_widths = c(1, 0.2), align = "h")
# Combine height for saving plot
combined_height <- max(rp_height, rs_height)
# Save the plot
setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/")
#ggsave("Supplementary_FigSXX_BubblePlot_byStage.tiff", 
       #combined.bubble.final, 
       #device = "tiff", width = 18.2, 
       #height = combined_height, 
       #units = "in", dpi = 300, 
       #compression = "lzw", bg = "white")










