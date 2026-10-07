############################################################################################################################################
#Title: "Rhizosphere Assembly: Sloan Neutral Model by Growth Stage and Stacked Bar Plot (Mean Relative Abundance) of Host-Selected ASVs" ###
#Author: "Ari Fina Bintarti"
#Date: "26-09-2026"
############################################################################################################################################
# Host-Selected ASVs = Overlapped/shared ASVs between Indicator Species analysis (IndVal.g) and Sloan Neutral Model = Above Prediction

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(ggnewscale)


#############################################################################
### Prepare Data Set for Plotting Stacked Bar Plots of Host-Selected Taxa ###
#############################################################################


# Read the updated Uniform Taxonomy SILVA 144
new.tax.unif <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/uniform_taxonomy_edit.csv", row.names = 1)
new.tax.unif <- rownames_to_column(new.tax.unif, var = "asvID")
# Read the updated Weighted Taxonomy SILVA 144 (Weighted = based on specific habitat/sampling location, here I used the weighted taxonomy of Plant-Rhizosphere Soil = https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/weighted/plant-rhizosphere)
new.tax.weight <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/weight_taxonomy_edit.csv", row.names = 1)
new.tax.weight <- rownames_to_column(new.tax.weight, var = "asvID")


# I. Rhizoplane

# Read the rhizoplane indicator stage ASVs data which also contain all information from the other non-indicator stages
rp.indic.above.all.ed <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizoplane_indic_above_all.csv", row.names = 1)
dim(rp.indic.above.all.ed) # 744
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
# change stage level
rp.indic.above.all.edit$Stage <- factor(
  rp.indic.above.all.edit$Stage, levels = c("Seedling","Early Vegetative","Late Vegetative",
                                          "Early Reproductive","Late Reproductive"))
# Make the Stage in two rows
rp.indic.above.all.edit$Stage <- factor(
  rp.indic.above.all.edit$Stage, levels = c("Seedling","Early Vegetative","Late Vegetative",
                                          "Early Reproductive","Late Reproductive"),
  labels = c("Seedling","Early\nVegetative","Late\nVegetative",
             "Early\nReproductive","Late\nReproductive"))
# Define stage order
stage_order <- c("Seedling","Early\nVegetative","Late\nVegetative","Early\nReproductive","Late\nReproductive")
# Add Compartment
rp.indic.above.all.edit <- rp.indic.above.all.edit %>%
  group_by(Weight.Family_asvID) %>%
  mutate(indicator_stage = Stage[stage_indicator][1]) %>%
  ungroup() %>%
  filter(stage_indicator == TRUE) %>%
  mutate(Compartment = "Rhizoplane")


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
# change stage level
rs.indic.above.all.edit$Stage <- factor(
  rs.indic.above.all.edit$Stage, levels = c("Seedling","Early Vegetative","Late Vegetative",
                                          "Early Reproductive","Late Reproductive"))
# Make the Stage in two rows
rs.indic.above.all.edit$Stage <- factor(rs.indic.above.all.edit$Stage,
                                      levels = c("Seedling","Early Vegetative","Late Vegetative","Early Reproductive","Late Reproductive"),
                                      labels = c("Seedling","Early\nVegetative","Late\nVegetative","Early\nReproductive","Late\nReproductive"))
# Define stage order
stage_order <- c("Seedling","Early\nVegetative","Late\nVegetative","Early\nReproductive","Late\nReproductive")
# Add Compartment
rs.indic.above.all.edit <- rs.indic.above.all.edit %>%
  group_by(Weight.Family_asvID) %>%
  mutate(indicator_stage = Stage[stage_indicator][1]) %>%
  ungroup() %>%
  filter(stage_indicator == TRUE) %>%
  mutate(Compartment = "Rhizosphere")

# Combine Rhizoplane and Rhizosphere data to make a common color palette
indicator_all <- bind_rows(rp.indic.above.all.edit, rs.indic.above.all.edit)
# Determine full color set for all selected ASVs in rhizoplane and rhizosphere (there are 24 Phyla in total)
all_phyla_bar <- sort(unique(c(indicator_all$Weight.Phylum)))
tableau32 <- c("#F28E2B","#FFCDA1","#324DA0","#C3DBFD","#59A14F","#C2EFB4","#E15759","#FF9D9A",
               "#BC6EB9","#D4A6C8","#7D4B4B","#D7B5A6","#D66982","#FABFD2","#79706E","#BAB0AC",
               "#B6992D","#C0B878","#0FCFC0","#99CFD1","#8942bd","#D7BFE3","#FDE333","#FEFDBE",
               "#4B0055","#BA4B8E","#00214E","#0072B4","#9BB306","#798233","#802A07","#A36B2B")
phylum_cols_bar <- setNames(tableau32, all_phyla_bar)
phylum_order_bar <- names(phylum_cols_bar)
# Calculate the abundance by compartment x indicator stage x phylum
indicator_phylum_abund <- indicator_all %>%
  group_by(Compartment, indicator_stage, Weight.Phylum) %>%
  summarise(mean_rel_abund = sum(mean_relabund_percent, na.rm = TRUE), .groups = "drop") %>%
  mutate(Compartment = factor(Compartment, levels = c("Rhizoplane", "Rhizosphere")),
         indicator_stage = factor(indicator_stage, levels = stage_order),
         Weight.Phylum = factor(Weight.Phylum, levels = phylum_order_bar))
# Calculate the number of the indicator ASVs
stage_counts <- indicator_all %>% 
  distinct(Compartment, indicator_stage, Weight.Family_asvID) %>%
  count(Compartment, indicator_stage, name = "n_ASV")
stage_totals <- indicator_phylum_abund %>%
  group_by(Compartment, indicator_stage) %>%
  summarise(total_abund = sum(mean_rel_abund), .groups = "drop") %>%
  left_join(stage_counts, by = c("Compartment", "indicator_stage"))

# Separate between Rhizoplane and Rhizosphere

# 1. Rhizoplane

# Keep only Rhizoplane
indicator_rhizoplane <- indicator_all %>%
  filter(Compartment == "Rhizoplane")
# Calculate the abundance by indicator stage x phylum
rp.indicator_phylum_abund <- indicator_rhizoplane %>%
  group_by(indicator_stage, Weight.Phylum) %>%
  summarise(mean_rel_abund = sum(mean_relabund_percent, na.rm = TRUE),
            .groups = "drop") %>%
  dplyr::mutate(indicator_stage = factor(indicator_stage, levels = stage_order),
                Weight.Phylum = factor(Weight.Phylum, levels = phylum_order_bar))
# Calculate the number of indicator ASVs
rp.stage_counts <- indicator_rhizoplane %>% 
  distinct(indicator_stage, Weight.Family_asvID) %>%
  count(indicator_stage, name = "n_ASV")
# Join with the total abundance
rp.stage_totals <- rp.indicator_phylum_abund %>%
  group_by(indicator_stage) %>%
  summarise(total_abund = sum(mean_rel_abund),
            .groups = "drop") %>%
  left_join(rp.stage_counts, by = "indicator_stage")
# Rhizoplane stacked bar plot
rp.indicator_abundance_plot <- ggplot(rp.indicator_phylum_abund,
                                      aes(x = indicator_stage, y = mean_rel_abund, fill = Weight.Phylum)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = phylum_cols_bar,
                    breaks = phylum_order_bar, drop = TRUE) +
  labs(title="A6. Host-Selected Taxa", x = NULL, y = "Mean Relative Abundance (%)", fill = "Phylum") +
  theme_bw() +
  theme(plot.title = element_text(size=15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(color = "transparent", size = 14, vjust = -0.75),
        axis.title.y = element_text(size = 13),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.width = unit(0, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 15, l = 5.5))+
  #panel.border = element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.5)) +
  geom_text(data = rp.stage_totals, aes(x = indicator_stage, y = total_abund, label = paste0("n = ", n_ASV)),
            inherit.aes = FALSE, vjust = -0.4, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 5.5)))


# 2. Rhizosphere

# Keep only Rhizosphere
indicator_rhizosphere <- indicator_all %>%
  filter(Compartment == "Rhizosphere")
# Calculate the abundance by indicator stage x phylum
rs.indicator_phylum_abund <- indicator_rhizosphere %>%
  group_by(indicator_stage, Weight.Phylum) %>%
  summarise(mean_rel_abund = sum(mean_relabund_percent, na.rm = TRUE),
            .groups = "drop") %>%
  dplyr::mutate(indicator_stage = factor(indicator_stage, levels = stage_order),
                Weight.Phylum = factor(Weight.Phylum, levels = phylum_order_bar))
# Calculate the number of indicator ASVs
rs.stage_counts <- indicator_rhizosphere %>% 
  distinct(indicator_stage, Weight.Family_asvID) %>%
  count(indicator_stage, name = "n_ASV")
# Join with the total abundance
rs.stage_totals <- rs.indicator_phylum_abund %>%
  group_by(indicator_stage) %>%
  summarise(total_abund = sum(mean_rel_abund),
            .groups = "drop") %>%
  left_join(rs.stage_counts, by = "indicator_stage")
# Rhizopshere stacked bar plot
rs.indicator_abundance_plot <- ggplot(rs.indicator_phylum_abund,
                                      aes(x = indicator_stage, y = mean_rel_abund, fill = Weight.Phylum)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = phylum_cols_bar,
                    breaks = phylum_order_bar, drop = TRUE) +
  labs(title="B6. Host-Selected Taxa", x = NULL, y = "Mean Relative Abundance (%)", fill = "Phylum") +
  theme_bw() +
  theme(plot.title = element_text(size=15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.width = unit(0, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.grid = element_blank(),
        #panel.border = element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.5),
        plot.margin = margin(t = 5.5, r = 5.5, b = 15, l = 5.5)) +
  geom_text(data = rs.stage_totals, aes(x = indicator_stage, y = total_abund, label = paste0("n = ", n_ASV)),
            inherit.aes = FALSE, vjust = -0.4, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 5.5)))



################################################################
### Prepare Data Set for the Sloan Neutral Model All Figures ###
################################################################

# I. Rhizoplane 

# 1. Seedling
# Load the Seedling Neutral Model and Indicator Species complete data set
rplane_stage1_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Seedling_complete_SNCM_Indic.csv", row.names = 1)
rplane_stage1_nm.rel.indic <- rownames_to_column(rplane_stage1_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rplane_stage1_nm.rel.indic <- rplane_stage1_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_stage1_nm.rel.indic <- rplane_stage1_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_stage1_nm.rel.indic <- rplane_stage1_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_stage1_nm.rel.indic.df <- rplane_stage1_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_stage1_nm.rel.indic.df$point_class <- ifelse(rplane_stage1_nm.rel.indic.df$above_CI == TRUE,
                                                    "Above",ifelse(rplane_stage1_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_stage1_nm.rel.indic.df$point_class <- factor(rplane_stage1_nm.rel.indic.df$point_class,
                                                    levels = c("Above", "Neutral", "Below"))

# 2.) Early Vegetative
# Load the Early Vegetative Neutral Model and Indicator Species complete data set
rplane_stage2_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Early_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
rplane_stage2_nm.rel.indic <- rownames_to_column(rplane_stage2_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rplane_stage2_nm.rel.indic <- rplane_stage2_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_stage2_nm.rel.indic <- rplane_stage2_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_stage2_nm.rel.indic <- rplane_stage2_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_stage2_nm.rel.indic.df <- rplane_stage2_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_stage2_nm.rel.indic.df$point_class <- ifelse(rplane_stage2_nm.rel.indic.df$above_CI == TRUE,
                                                    "Above",ifelse(rplane_stage2_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_stage2_nm.rel.indic.df$point_class <- factor(rplane_stage2_nm.rel.indic.df$point_class,
                                                    levels = c("Above", "Neutral", "Below"))

# 3.) Late Vegetative
# Load the Late Vegetative Neutral Model and Indicator Species complete data set
rplane_stage3_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Late_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
rplane_stage3_nm.rel.indic <- rownames_to_column(rplane_stage3_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rplane_stage3_nm.rel.indic <- rplane_stage3_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_stage3_nm.rel.indic <- rplane_stage3_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_stage3_nm.rel.indic <- rplane_stage3_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_stage3_nm.rel.indic.df <- rplane_stage3_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_stage3_nm.rel.indic.df$point_class <- ifelse(rplane_stage3_nm.rel.indic.df$above_CI == TRUE,
                                                    "Above",ifelse(rplane_stage3_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_stage3_nm.rel.indic.df$point_class <- factor(rplane_stage3_nm.rel.indic.df$point_class,
                                                    levels = c("Above", "Neutral", "Below"))

# 4.) Early Reproductive
# Load the Early Reproductive Neutral Model and Indicator Species complete data set
rplane_stage4_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Early_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
rplane_stage4_nm.rel.indic <- rownames_to_column(rplane_stage4_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rplane_stage4_nm.rel.indic <- rplane_stage4_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_stage4_nm.rel.indic <- rplane_stage4_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_stage4_nm.rel.indic <- rplane_stage4_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_stage4_nm.rel.indic.df <- rplane_stage4_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_stage4_nm.rel.indic.df$point_class <- ifelse(rplane_stage4_nm.rel.indic.df$above_CI == TRUE,
                                                    "Above",ifelse(rplane_stage4_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_stage4_nm.rel.indic.df$point_class <- factor(rplane_stage4_nm.rel.indic.df$point_class,
                                                    levels = c("Above", "Neutral", "Below"))

# 5.) Late Reproductive
# Load the Late Reproductive Neutral Model and Indicator Species complete data set
rplane_stage5_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Late_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
rplane_stage5_nm.rel.indic <- rownames_to_column(rplane_stage5_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rplane_stage5_nm.rel.indic <- rplane_stage5_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_stage5_nm.rel.indic <- rplane_stage5_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_stage5_nm.rel.indic <- rplane_stage5_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_stage5_nm.rel.indic.df <- rplane_stage5_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_stage5_nm.rel.indic.df$point_class <- ifelse(rplane_stage5_nm.rel.indic.df$above_CI == TRUE,
                                                    "Above",ifelse(rplane_stage5_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_stage5_nm.rel.indic.df$point_class <- factor(rplane_stage5_nm.rel.indic.df$point_class,
                                                    levels = c("Above", "Neutral", "Below"))


# II. Rhizosphere 

# 1. Seedling
# Load the Seedling Neutral Model and Indicator Species complete data set
rsphere_stage1_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Seedling_complete_SNCM_Indic.csv", row.names = 1)
rsphere_stage1_nm.rel.indic <- rownames_to_column(rsphere_stage1_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rsphere_stage1_nm.rel.indic <- rsphere_stage1_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_stage1_nm.rel.indic <- rsphere_stage1_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_stage1_nm.rel.indic <- rsphere_stage1_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_stage1_nm.rel.indic.df <- rsphere_stage1_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_stage1_nm.rel.indic.df$point_class <- ifelse(rsphere_stage1_nm.rel.indic.df$above_CI == TRUE,
                                                     "Above",ifelse(rsphere_stage1_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_stage1_nm.rel.indic.df$point_class <- factor(rsphere_stage1_nm.rel.indic.df$point_class,
                                                     levels = c("Above", "Neutral", "Below"))

# 2.) Early Vegetative
# Load the Early Vegetative Neutral Model and Indicator Species complete data set
rsphere_stage2_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Early_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
rsphere_stage2_nm.rel.indic <- rownames_to_column(rsphere_stage2_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rsphere_stage2_nm.rel.indic <- rsphere_stage2_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_stage2_nm.rel.indic <- rsphere_stage2_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_stage2_nm.rel.indic <- rsphere_stage2_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_stage2_nm.rel.indic.df <- rsphere_stage2_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_stage2_nm.rel.indic.df$point_class <- ifelse(rsphere_stage2_nm.rel.indic.df$above_CI == TRUE,
                                                     "Above",ifelse(rsphere_stage2_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_stage2_nm.rel.indic.df$point_class <- factor(rsphere_stage2_nm.rel.indic.df$point_class,
                                                     levels = c("Above", "Neutral", "Below"))

# 3.) Late Vegetative
# Load the Late Vegetative Neutral Model and Indicator Species complete data set
rsphere_stage3_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Late_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
rsphere_stage3_nm.rel.indic <- rownames_to_column(rsphere_stage3_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rsphere_stage3_nm.rel.indic <- rsphere_stage3_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_stage3_nm.rel.indic <- rsphere_stage3_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_stage3_nm.rel.indic <- rsphere_stage3_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_stage3_nm.rel.indic.df <- rsphere_stage3_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_stage3_nm.rel.indic.df$point_class <- ifelse(rsphere_stage3_nm.rel.indic.df$above_CI == TRUE,
                                                     "Above",ifelse(rsphere_stage3_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_stage3_nm.rel.indic.df$point_class <- factor(rsphere_stage3_nm.rel.indic.df$point_class,
                                                     levels = c("Above", "Neutral", "Below"))

# 4.) Early Reproductive
# Load the Early Reproductive Neutral Model and Indicator Species complete data set
rsphere_stage4_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Early_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
rsphere_stage4_nm.rel.indic <- rownames_to_column(rsphere_stage4_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rsphere_stage4_nm.rel.indic <- rsphere_stage4_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_stage4_nm.rel.indic <- rsphere_stage4_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_stage4_nm.rel.indic <- rsphere_stage4_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_stage4_nm.rel.indic.df <- rsphere_stage4_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_stage4_nm.rel.indic.df$point_class <- ifelse(rsphere_stage4_nm.rel.indic.df$above_CI == TRUE,
                                                     "Above",ifelse(rsphere_stage4_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_stage4_nm.rel.indic.df$point_class <- factor(rsphere_stage4_nm.rel.indic.df$point_class,
                                                     levels = c("Above", "Neutral", "Below"))

# 5.) Late Reproductive
# Load the Late Reproductive Neutral Model and Indicator Species data set
rsphere_stage5_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Late_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
rsphere_stage5_nm.rel.indic <- rownames_to_column(rsphere_stage5_nm.rel.indic, var = "asvID")
# Re-name the old taxonomy SILVA 138 column name
rsphere_stage5_nm.rel.indic <- rsphere_stage5_nm.rel.indic %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_stage5_nm.rel.indic <- rsphere_stage5_nm.rel.indic %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_stage5_nm.rel.indic <- rsphere_stage5_nm.rel.indic %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_stage5_nm.rel.indic.df <- rsphere_stage5_nm.rel.indic %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_stage5_nm.rel.indic.df$point_class <- ifelse(rsphere_stage5_nm.rel.indic.df$above_CI == TRUE,
                                                     "Above",ifelse(rsphere_stage5_nm.rel.indic.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_stage5_nm.rel.indic.df$point_class <- factor(rsphere_stage5_nm.rel.indic.df$point_class,
                                                     levels = c("Above", "Neutral", "Below"))

### Make Color Palette for The Host-Selected Taxa (These taxa were selected based on the overlapping/shared taxa between Sloan Neutral Model = Above Prediction and Indicator Species analysis)

# Load shared taxa of Sloan Neutral-selected and species indicator taxa in all growth stages
# 1.Rhizoplane:
rp.above_Ind.edtax <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Rhizoplane_Neutral_Indicator_SharedASVs.csv", row.names = 1)
# Check how many Phylum and what are they
unique(na.omit(rp.above_Ind.edtax$Unif.Phylum)) # 28 Phyla
unique(na.omit(rp.above_Ind.edtax$Weight.Phylum)) # 27 Phyla

# 2.Rhizosphere:
rsphere.above_Ind.edtax <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Rhizosphere_Neutral_Indicator_SharedASVs.csv", row.names = 1)
# Check how many Phylum and what are they
unique(na.omit(rsphere.above_Ind.edtax$Unif.Phylum)) # 24 Phyla
unique(na.omit(rsphere.above_Ind.edtax$Weight.Phylum)) # 24 Phyla

# In this analysis, we will use the weighted taxonomy SILVA 144
# Determine color for rhizoplane and rhizosphere Phyla 
all_phyla <- sort(unique(c(rp.above_Ind.edtax$Weight.Phylum, 
                           rsphere.above_Ind.edtax$Weight.Phylum))) # 32 total phyla in all compartments
# Common color palette for all phyla
tableau32 <- c("#F28E2B","#FFCDA1","#324DA0","#C3DBFD","#59A14F","#C2EFB4","#E15759","#FF9D9A",
               "#BC6EB9","#D4A6C8","#7D4B4B","#D7B5A6","#D66982","#FABFD2","#79706E","#BAB0AC",
               "#B6992D","#C0B878","#0FCFC0","#99CFD1","#8942bd","#D7BFE3","#FDE333","#FEFDBE",
               "#4B0055","#BA4B8E","#00214E","#0072B4","#9BB306","#798233","#802A07","#A36B2B")
phylum_cols <- setNames(tableau32, all_phyla)
# save the Phylum color palette
#saveRDS(phylum_cols, file = "/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/phylum_cols.rds")

### Plotting Sloan Neutral Model by growth stage and indicate the host-selected taxa

# I. Rhizoplane

# 1.) Seedling
# Make the SNCM Plot in Seedling
rp.seed.SNCM.indic.plot <- ggplot(rplane_stage1_nm.rel.indic.df, 
                                  aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" ="#44AA99"), 
                     labels = c("Above" = "Above (A)",
                                "Neutral" = "Neutral (N)",
                                "Below" = "Below (B)"),
                     guide = guide_legend(override.aes = list(shape = 16, fill =NA, size=3))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rplane_stage1_nm.rel.indic.df, stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="A1. Seedling",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y =element_text(size=13),
        axis.text=element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        #legend.position = c(0.19, 0.72),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.71\nm = 1.5e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 30\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Seedling
rp.stage1.predfit_count <- rplane_stage1_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(2.4%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(96.1%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.6%)")))
#label = paste0(predictions.fit_class, "\n", n, " (", round(percentage, 1), "%)"))
rp.stage1.predfit_count$predictions.label <- factor(rp.stage1.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rp.stage1.predfit_count$pos = (cumsum(c(0, rp.stage1.predfit_count$n)) + 
                                 c(rp.stage1.predfit_count$n / 2, .01))[1:nrow(rp.stage1.predfit_count)]
# Make Neutral Model prediction class pie chart in Seedling
rp.seed.SNCM.pie <- ggplot(rp.stage1.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rp.stage1.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rp.stage1.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) +
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rp.seed.SNCM.indic.pie.plot <- ggdraw(rp.seed.SNCM.indic.plot) +
  draw_plot(rp.seed.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)


# 2.) Early Vegetative
# Make the SNCM Plot in Early Vegetative
rp.early.veget.SNCM.indic.plot <- ggplot(rplane_stage2_nm.rel.indic.df, 
                                         aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rplane_stage2_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="A2. Early Vegetative",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text = element_text(size=12),
        #axis.title.y = element_text(color = "transparent", size = 14, vjust = -0.75),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.35, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.75\nm = 2e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  # annotate("text",x = -7,y = 0.75, label = "Indicator = 25\nTaxa",
  # size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in early Vegetative
rp.stage2.predfit_count <- rplane_stage2_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.9%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(96.5%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.6%)")))
#label = paste0(predictions.fit_class, "\n", n, " (", round(percentage, 1), "%)"))
rp.stage2.predfit_count$predictions.label <- factor(rp.stage2.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rp.stage2.predfit_count$pos = (cumsum(c(0, rp.stage2.predfit_count$n)) + 
                                 c(rp.stage2.predfit_count$n / 2, .01))[1:nrow(rp.stage2.predfit_count)]
# Make Neutral Model prediction class pie chart in Early Vegetative
rp.early.veget.SNCM.pie <- ggplot(rp.stage2.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rp.stage2.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rp.stage2.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) +
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rp.early.veget.SNCM.indic.pie.plot <- ggdraw(rp.early.veget.SNCM.indic.plot) +
  draw_plot(rp.early.veget.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3) 

# 3.) Late Vegetative
# Make the SNCM Plot in Late Vegetative
rp.late.veget.SNCM.indic.plot <- ggplot(rplane_stage3_nm.rel.indic.df, 
                                        aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rplane_stage3_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="A3. Late Vegetative",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        #axis.title.x = element_text(color = "transparent", size = 20, vjust = -0.75),
        axis.title.x = element_blank(),
        axis.title.y =element_text(size=13),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.37, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.76\nm = 2e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 54\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Late Vegetative
rp.stage3.predfit_count <- rplane_stage3_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.5%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(97.2%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.3%)")))
rp.stage3.predfit_count$predictions.label <- factor(rp.stage3.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rp.stage3.predfit_count$pos = (cumsum(c(0, rp.stage3.predfit_count$n)) + 
                                 c(rp.stage3.predfit_count$n / 2, .01))[1:nrow(rp.stage3.predfit_count)]
# Make Neutral Model prediction class pie chart in Late Vegetative
rp.late.veget.SNCM.pie <- ggplot(rp.stage3.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rp.stage3.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rp.stage3.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4,lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rp.late.veget.SNCM.indic.pie.plot <- ggdraw(rp.late.veget.SNCM.indic.plot) +
  draw_plot(rp.late.veget.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3) +
  draw_label("log10(mean relative abundance)",x = 0.5,y = 0.03,size = 13)

# 4.) Early Reproductive
# Make the SNCM Plot in Early Reproductive
rp.early.reprod.SNCM.indic.plot <- ggplot(rplane_stage4_nm.rel.indic.df, 
                                          aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rplane_stage4_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="A4. Early Reproductive",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "transparent", size = 13, vjust = -0.75),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.37, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.77\nm = 1.8e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 26\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Early Reproductive
rp.stage4.predfit_count <- rplane_stage4_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.8%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(96.8%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.3%)")))
rp.stage4.predfit_count$predictions.label <- factor(rp.stage4.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rp.stage4.predfit_count$pos = (cumsum(c(0, rp.stage4.predfit_count$n)) + 
                                 c(rp.stage4.predfit_count$n / 2, .01))[1:nrow(rp.stage4.predfit_count)]
# Make Neutral Model prediction class pie chart in Early Reproductive
rp.early.reprod.SNCM.pie <- ggplot(rp.stage4.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rp.stage4.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rp.stage4.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rp.early.reprod.SNCM.indic.pie.plot <- ggdraw(rp.early.reprod.SNCM.indic.plot) +
  draw_plot(rp.early.reprod.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)

# 5.) Late Reproductive
# Make the SNCM Plot in Late Reproductive
rp.late.reprod.SNCM.indic.plot <- ggplot(rplane_stage5_nm.rel.indic.df, 
                                         aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     labels = c("Above" = "Above pred. (A)",
                                "Neutral" = "Neutral (N)",
                                "Below" = "Below pred. (B)"),
                     guide = guide_legend(override.aes = list(shape = 16, fill =NA, size=3, stroke =1))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rplane_stage5_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3.2) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="A5. Late Reproductive",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(size=20, vjust = -0.75),
        axis.title.y = element_text(color = "transparent", size = 13, vjust = -0.75),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=9.8),
        #legend.position = c(0.85, 0.12),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.width = unit(0, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.68\nm = 1.4e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 24\nTaxa",
  # size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")

# Create count data of the Neutral Model prediction class in Late Reproductive
rp.stage5.predfit_count <- rplane_stage5_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.8%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(96.8%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.4%)")))
rp.stage5.predfit_count$predictions.label <- factor(rp.stage5.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rp.stage5.predfit_count$pos = (cumsum(c(0, rp.stage5.predfit_count$n)) + 
                                 c(rp.stage5.predfit_count$n / 2, .01))[1:nrow(rp.stage5.predfit_count)]
# Make Neutral Model prediction class pie chart in Late Reproductive
rp.late.reprod.SNCM.pie <- ggplot(rp.stage5.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rp.stage5.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rp.stage5.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rp.late.reprod.SNCM.indic.pie.plot <- ggdraw(rp.late.reprod.SNCM.indic.plot) +
  draw_plot(rp.late.reprod.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)+
  draw_label("log10(mean relative abundance)",x = 0.55,y = 0.03,size = 13)

#ADDING STACKED BARPLOT HERE!!!
rp.indicator_abundance_plot_ggdr <- ggdraw(rp.indicator_abundance_plot) 

# Combine rhizoplane 
rp.all.SNCM.indic.pie.plot <- cowplot::plot_grid(rp.seed.SNCM.indic.pie.plot,
                                                 rp.early.reprod.SNCM.indic.pie.plot,
                                                 rp.early.veget.SNCM.indic.pie.plot,
                                                 rp.late.reprod.SNCM.indic.pie.plot,
                                                 rp.late.veget.SNCM.indic.pie.plot,
                                                 rp.indicator_abundance_plot_ggdr,
                                                 ncol = 2, align = "v",
                                                 axis = "lr")
# Add the common title of Rhizoplane
rp.all.SNCM.indic.pie.plot_w_title <- ggdraw() + 
  draw_plot(plot_grid(ggdraw() + draw_label("A. Rhizoplane", x = 0, y = 0.95, size = 23, hjust=0, vjust=1),
  rp.all.SNCM.indic.pie.plot,ncol = 1,rel_heights = c(0.05, 1)), x=0, y=0.01, width = 0.99,height = 0.99)


# II. Rhizosphere

# 1.) Seedling
# Make the SNCM Plot in Seedling
rs.seed.SNCM.indic.plot <- ggplot(rsphere_stage1_nm.rel.indic.df, 
                                  aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rsphere_stage1_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="B1. Seedling",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.37, 0.89),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.85\nm = 1.8e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 29\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Seedling
rs.stage1.predfit_count <- rsphere_stage1_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.5%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(97.1%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.4%)")))
#label = paste0(predictions.fit_class, "\n", n, " (", round(percentage, 1), "%)"))
rs.stage1.predfit_count$predictions.label <- factor(rs.stage1.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rs.stage1.predfit_count$pos = (cumsum(c(0, rs.stage1.predfit_count$n)) + 
                                 c(rs.stage1.predfit_count$n / 2, .01))[1:nrow(rs.stage1.predfit_count)]
# Make Neutral Model prediction class pie chart in Seedling
rs.seed.SNCM.pie <- ggplot(rs.stage1.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rs.stage1.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rs.stage1.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) +
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rs.seed.SNCM.indic.pie.plot <- ggdraw(rs.seed.SNCM.indic.plot) +
  draw_plot(rs.seed.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)

# 2.) Early Vegetative
# Make the SNCM Plot in Early Vegetative
rs.early.veget.SNCM.indic.plot <- ggplot(rsphere_stage2_nm.rel.indic.df, 
                                         aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rsphere_stage2_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="B2. Early Vegetative",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text = element_text(size=12),
        #axis.title.y = element_text(color = "transparent", size = 14, vjust = -0.75),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.35, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.79\nm = 1.7e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 32\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in early Vegetative
rs.stage2.predfit_count <- rsphere_stage2_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.5%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(97.1%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.4%)")))
#label = paste0(predictions.fit_class, "\n", n, " (", round(percentage, 1), "%)"))
rs.stage2.predfit_count$predictions.label <- factor(rs.stage2.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rs.stage2.predfit_count$pos = (cumsum(c(0, rs.stage2.predfit_count$n)) + 
                                 c(rs.stage2.predfit_count$n / 2, .01))[1:nrow(rs.stage2.predfit_count)]
# Make Neutral Model prediction class pie chart in Early Vegetative
rs.early.veget.SNCM.pie <- ggplot(rs.stage2.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rs.stage2.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rs.stage2.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) +
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rs.early.veget.SNCM.indic.pie.plot <- ggdraw(rs.early.veget.SNCM.indic.plot) +
  draw_plot(rs.early.veget.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)

# 3.) Late Vegetative
# Make the SNCM Plot in Late Vegetative
rs.late.veget.SNCM.indic.plot <- ggplot(rsphere_stage3_nm.rel.indic.df, 
                                        aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rsphere_stage3_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="B3. Late Vegetative",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        #axis.title.x = element_text(color = "transparent", size = 20, vjust = -0.75),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=13),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.37, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.84\nm = 1.8e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 18\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Late Vegetative
rs.stage3.predfit_count <- rsphere_stage3_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.7%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(96.9%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.4%)")))
rs.stage3.predfit_count$predictions.label <- factor(rs.stage3.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rs.stage3.predfit_count$pos = (cumsum(c(0, rs.stage3.predfit_count$n)) + 
                                 c(rs.stage3.predfit_count$n / 2, .01))[1:nrow(rs.stage3.predfit_count)]
# Make Neutral Model prediction class pie chart in Late Vegetative
rs.late.veget.SNCM.pie <- ggplot(rs.stage3.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rs.stage3.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rs.stage3.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4,lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rs.late.veget.SNCM.indic.pie.plot <- ggdraw(rs.late.veget.SNCM.indic.plot) +
  draw_plot(rs.late.veget.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3) +
  draw_label("log10(mean relative abundance)",x = 0.56,y = 0.03,size = 13)

# 4.) Early Reproductive
# Make the SNCM Plot in Early Reproductive
rs.early.reprod.SNCM.indic.plot <- ggplot(rsphere_stage4_nm.rel.indic.df, 
                                          aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rsphere_stage4_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="B4. Early Reproductive",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12),
        axis.title.y = element_text(color = "transparent", size = 13, vjust = -0.75),
        #axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = "none",
        #legend.position = c(0.37, 0.89),
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.82\nm = 1.6e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 14\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")
# Create count data of the Neutral Model prediction class in Early Reproductive
rs.stage4.predfit_count <- rsphere_stage4_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.6%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(97.2%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.3%)")))
rs.stage4.predfit_count$predictions.label <- factor(rs.stage4.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rs.stage4.predfit_count$pos = (cumsum(c(0, rs.stage4.predfit_count$n)) + 
                                 c(rs.stage4.predfit_count$n / 2, .01))[1:nrow(rs.stage4.predfit_count)]
# Make Neutral Model prediction class pie chart in Early Reproductive
rs.early.reprod.SNCM.pie <- ggplot(rs.stage4.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rs.stage4.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rs.stage4.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rs.early.reprod.SNCM.indic.pie.plot <- ggdraw(rs.early.reprod.SNCM.indic.plot) +
  draw_plot(rs.early.reprod.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3) 
#draw_label("log10(mean relative abundance)",x = 0.56,y = 0,size = 14)

# 5.) Late Reproductive
# Make the SNCM Plot in Late Reproductive
rs.late.reprod.SNCM.indic.plot <- ggplot(rsphere_stage5_nm.rel.indic.df, 
                                         aes(x = log_mean, y = occupancy, color = point_class)) +
  geom_point(shape = 21,fill = "white",alpha = 0.7,size = 2) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     guide = guide_legend(override.aes = list(shape = 16,fill = NA, size=4))) +
  # Reset color scale
  ggnewscale::new_scale_fill() +
  # Highlight indicator taxa and above prediction, colored by taxonomy
  geom_point(data = subset(rsphere_stage5_nm.rel.indic.df,stage_indicator == TRUE & above_CI == TRUE),
             aes(fill = Weight.Phylum),shape = 21,color = "black",size = 3) +
  scale_fill_manual(values = phylum_cols) +
  geom_line(color = "black",linewidth = 0.5, 
            aes(y = predictions.freq.pred,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.upr,x = log10(predictions.p)),alpha = 0.5) +
  geom_line(color = "black",linetype = "dashed",linewidth = 0.5,
            aes(y = predictions.pred.lwr,x = log10(predictions.p)),alpha = 0.5)+
  labs(title="B5. Late Reproductive",x="log10(mean relative abundance)", y="Occupancy", fill="Phylum") +
  theme_bw()+
  theme(plot.title = element_text(size=15),
        axis.title.x = element_blank(),
        #axis.title.x = element_text(size=20, vjust = -0.75),
        axis.title.y = element_text(color = "transparent", size = 13, vjust = -0.75),
        axis.text = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        #legend.position = c(0.37, 0.89),
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        panel.background = element_rect(fill = "white", colour = "white"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 15, 5.5))+
  annotate("text",x = -7,y = 0.94, label = "R² = 0.82\nm = 1.7e-4",
           size = 4, hjust=0, fontface='italic', lineheight=0.8)+
  #annotate("text",x = -7,y = 0.75, label = "Indicator = 38\nTaxa",
  #size = 4, hjust=0, fontface='italic', lineheight=0.8) +
  guides(fill = "none")

# Create count data of the Neutral Model prediction class in Late Reproductive
rs.stage5.predfit_count <- rsphere_stage5_nm.rel.indic.df %>%
  count(predictions.fit_class) %>%
  mutate(percentage = n / sum(n) * 100,
         label = case_when(predictions.fit_class == "Above prediction" ~ paste0("A:", n, "\n(1.5%)"),
                           predictions.fit_class == "As predicted" ~ paste0("N:", n, "\n(97.2%)"),
                           predictions.fit_class == "Below prediction" ~ paste0("B:", n, "\n(1.4%)")))
rs.stage5.predfit_count$predictions.label <- factor(rs.stage5.predfit_count$predictions.fit_class, levels = c("Above prediction", "As predicted", "Below prediction"),
                                                    labels = c("Above", "Neutral", "Below"))
rs.stage5.predfit_count$pos = (cumsum(c(0, rs.stage5.predfit_count$n)) + 
                                 c(rs.stage5.predfit_count$n / 2, .01))[1:nrow(rs.stage5.predfit_count)]
# Make Neutral Model prediction class pie chart in Late Reproductive
rs.late.reprod.SNCM.pie <- ggplot(rs.stage5.predfit_count, aes(x = "", y = n, fill = predictions.label)) +
  geom_col(position = position_stack(reverse = TRUE), 
           show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Above" = "#DDCC77",
    "Neutral" = "lightgray",
    "Below" = "#44AA99"))+
  theme_void()+
  theme(legend.position="none")+
  geom_text_repel(data = subset(rs.stage5.predfit_count, predictions.label != "Neutral"),
                  aes(x = 1.4, y = pos, label = label, color=predictions.label), 
                  size=4, nudge_x = .3, 
                  segment.size = .5, 
                  show.legend = FALSE, 
                  lineheight=0.8)+
  geom_text(data = subset(rs.stage5.predfit_count, predictions.label == "Neutral"),
            aes(x = 1,y = pos,label = label,color = predictions.label),
            size = 4, lineheight=0.8,
            show.legend = FALSE) +
  scale_color_manual(values = c(
    "Above" = "#BFA94F",
    "Neutral" = "gray40",
    "Below" = "#44AA99")) + 
  coord_polar('y') +
  theme_void()
# Overlay pie chart inside the SNCM Plot 
rs.late.reprod.SNCM.indic.pie.plot <- ggdraw(rs.late.reprod.SNCM.indic.plot) +
  draw_plot(rs.late.reprod.SNCM.pie,x = 0.15,y = 0.25,width = 0.3,height = 0.3)+
  draw_label("log10(mean relative abundance)",x = 0.55,y = 0.03,size = 13)

#ADDING STACKED BARPLOT HERE!!!
rs.indicator_abundance_plot_ggdr <- ggdraw(rs.indicator_abundance_plot) 

# Combine all rhizosphere plots with the stacked barplot 
rs.all.SNCM.indic.pie.plot <- cowplot::plot_grid(rs.seed.SNCM.indic.pie.plot,
                                                 rs.early.reprod.SNCM.indic.pie.plot,
                                                 rs.early.veget.SNCM.indic.pie.plot,
                                                 rs.late.reprod.SNCM.indic.pie.plot,
                                                 rs.late.veget.SNCM.indic.pie.plot,
                                                 rs.indicator_abundance_plot_ggdr,
                                                 ncol = 2, align = "v",
                                                 axis = "lr")

# Add the common title of Rhizosphere
rs.all.SNCM.indic.pie.plot_w_title <- ggdraw()+
  draw_plot(plot_grid(ggdraw() + draw_label("B. Rhizosphere", x = 0, y = 0.95, size = 23, hjust=0, vjust=1),
                      rs.all.SNCM.indic.pie.plot,ncol = 1,rel_heights = c(0.05, 1)), x=0, y=0.01, width = 0.99,height = 0.99)


##############################################
### Combine All Rhizosphere and Rhizoplane ###
##############################################

# Make Dummy legend-only plot for Predictions
prediction_legend_dummy <- ggplot(data.frame(
  point_class = factor(c("Above", "Neutral", "Below"),levels = c("Above", "Neutral", "Below")),
  x = 1:3, y = 1:3),
  aes(x = x, y = y, color = point_class)) +
  geom_point(size = 4) +
  scale_color_manual(name = "Prediction",
                     values = c("Above" = "#DDCC77","Neutral" = "gray","Below" = "#44AA99"),
                     labels = c("Above" = "Above prediction (A)",
                                "Neutral" = "Neutral (N)",
                                "Below" = "Below prediction (B)"),
                     breaks = c("Above", "Neutral", "Below"),
                     guide = guide_legend(override.aes = list(shape = 21,fill = NA,size = 3, stroke=1.5))) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.height = unit(0.8, "cm"))

# Make Dummy legend-only plot for Phylum
dummy_plot <-  # Stacked bar plot
  indicator_abundance_plot <- ggplot(indicator_phylum_abund,
                                     aes(x = indicator_stage, y = mean_rel_abund, fill = Weight.Phylum)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Compartment, ncol = 1, strip.position = "top") +
  scale_fill_manual(values = phylum_cols_bar, breaks = phylum_order_bar, drop = TRUE) +
  labs(x = NULL, y = "Mean Relative Abundance (%)", fill = "Phylum") +
  #scale_y_continuous(expand = expansion(mult = c(, 0.09))) +
  theme_bw() +
  theme(axis.text.y = element_text(size=17),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=17),
        axis.title.x = element_blank(),
        axis.title.y =element_text(size=22),
        legend.title.position = "top",
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        strip.text.x = element_text(size = 25),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.5)) +
  geom_text(data = stage_totals, aes(x = indicator_stage, y = total_abund, label = paste0("n = ", n_ASV)),
            inherit.aes = FALSE, vjust = -0.4, size = 6) +
  guides(fill = guide_legend(ncol = 8, override.aes = list(size = 5.5)))


# Extract the legend as a grob (graphical object)
# Prediction legend
prediction_legend <- cowplot::get_legend(prediction_legend_dummy)
# Phylum legend
legend_grob <- cowplot::get_legend(dummy_plot)
# Create a centered canvas for the two legends
prediction_centered <- cowplot::plot_grid(NULL, prediction_legend, 
                                          NULL, ncol = 3, rel_widths = c(0.2, 1, 0.8))
phylum_centered <- cowplot::plot_grid(NULL, legend_grob, 
                                      NULL, ncol = 3, rel_widths = c(0.00, 1, 0.10))
# Put the two centered legends side-by-side
legends_combined <- cowplot::plot_grid(prediction_centered, phylum_centered,
                                       ncol = 2, rel_widths = c(0.2, 0.8))
# Combine the Rhizoplane and Rhizosphere Plots
SNCM.Indic.plot <- wrap_plots(rp.all.SNCM.indic.pie.plot_w_title, 
                              rs.all.SNCM.indic.pie.plot_w_title, ncol = 2)
# Final figure
SNCM.Indic.final.plot <- cowplot::plot_grid(SNCM.Indic.plot, legends_combined, 
                                            ncol = 1, rel_heights = c(1, 0.17), align = "v")
# Save the plot
setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/")
#ggsave("FigXX_SloanNeutralModel_byStage.tiff",
       #SNCM.Indic.final.plot, device = "tiff",
       #width = 18.5, height =13, 
       #units= "in", dpi = 300,
       #compression="lzw", bg= "white")























