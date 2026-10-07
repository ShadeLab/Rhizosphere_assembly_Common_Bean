############################################################################################################################################
#title: "Rhizosphere Assembly: Indicator Species and Sloan Neutral Community Model Analyses"
#author: "Ari Fina Bintarti"
#date: "30-09-2026"
############################################################################################################################################

library(phyloseq)
library(readxl)
library(stringr)
library(writexl)
library(remotes)
library(minpack.lm)
library(Hmisc)
library(stats4)
#remotes::install_github("DanielSprockett/reltools", force = T)
library(reltools)
#install.packages("indicspecies")
library(indicspecies)
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#remotes::install_github("david-barnett/microViz")
library(microViz)
#install.packages("BRCore")
library(BRCore)
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)


# Load metadata with "stage aggregate"
setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/")
metadata2 <- read.csv("metadata2.csv", header=TRUE)
metadata2 <- column_to_rownames(metadata2, var ="SampleID")
# Load the multi-rarefied phyloseq object from Sam Barnett
RA_phyloseq_multrare <- readRDS("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/RA_phyloseq_multirarefied_Sam_16S.rds")
RA_phyloseq_multrare # 146,223 taxa and 269 samples
# The "RA_phyloseq_multrare" is a relative abundance ASV table (not raw count, not integer) as shown by the colSums:
sort(colSums(otu_table(RA_phyloseq_multrare), na.rm = FALSE, dims = 1), decreasing = T)
# Multiply everything in it by 20,000. This does not change any of the actual relative abundances but makes everything an integer. 
# This number is based on the fact that the OTU table calculated relative abundance across 20,000 rarefied reads and 100 replicate rarefactions. 
set.seed(13)
otu_table(RA_phyloseq_multrare) <- otu_table(RA_phyloseq_multrare)*20000*100 # 100 = the iteration number
sort(colSums(otu_table(RA_phyloseq_multrare), na.rm = FALSE, dims = 1), decreasing = T)
# Build a new phyloseq object using metadata2 that has "stage_aggregate"
asv <- otu_table(RA_phyloseq_multrare)
meta <- sample_data(metadata2)
tax <- tax_table(RA_phyloseq_multrare)
# Make a new phyloseq object
RAphyloseq <- merge_phyloseq(meta, tax, asv)
RAphyloseq # 146,223 taxa and 269 samples
any(taxa_sums(RAphyloseq) == 0) #TRUE
# Remove ASVs that are absent from all samples (i.e. all zeros), and keep ASVs even if they occur in only one sample.
RAps_nonzero <- prune_taxa(taxa_sums(RAphyloseq) > 0, RAphyloseq)
RAps_nonzero # 143,715 taxa and 269 samples
any(taxa_sums(RAps_nonzero) == 0) #FALSE
# Make data frames from the new phyloseq object
RAsampledata <- data.frame(sample_data(RAps_nonzero))
RA_ASVtable <- data.frame(otu_table(RAps_nonzero))
RA_taxtable <- data.frame(tax_table(RAps_nonzero))
# Subset by compartments: Rhizosphere and Rhizoplane
rsphere_ps <- RAps_nonzero %>% ps_filter(compartment == "rhizosphere") #85345 taxa and 136 samples
any(taxa_sums(rsphere_ps) == 0) #FALSE
rplane_ps <- RAps_nonzero %>% ps_filter(compartment == "rhizoplane") #65206 taxa and 133 samples
any(taxa_sums(rplane_ps) == 0) #FALSE
# Change integer and character into factor
metadata2[sapply(metadata2, is.integer)] <- lapply(metadata2[sapply(metadata2, is.integer)], as.factor)
metadata2[sapply(metadata2, is.character)] <- lapply(metadata2[sapply(metadata2, is.character)], as.factor)
metadata2$growthstage_aggregate_categorical
# Subset the phyloseq object by growth stage aggregate categorical
# 1.) Rhizosphere: seedling
rs_seedl_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "seedling")
any(taxa_sums(rs_seedl_ps) == 0)
# 2.) Rhizosphere: early vegetative
rs_earveg_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "early_vegetative")
any(taxa_sums(rs_earveg_ps) == 0)
# 3.) Rhizosphere: late vegetative
rs_lateveg_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "late_vegetative")
any(taxa_sums(rs_lateveg_ps) == 0)
# 4.) Rhizosphere: early reproductive
rs_earrep_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "early_reproductive")
any(taxa_sums(rs_earrep_ps) == 0)
# 5.) Rhizosphere: late reproductive
rs_laterep_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "late_reproductive")
any(taxa_sums(rs_laterep_ps) == 0)
# 1.) Rhizoplane: seedling
rp_seedl_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "seedling")
any(taxa_sums(rp_seedl_ps) == 0)
# 2.) Rhizoplane: early vegetative
rp_earveg_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "early_vegetative")
any(taxa_sums(rp_earveg_ps) == 0)
# 3.) Rhizoplane: late vegetative
rp_lateveg_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "late_vegetative")
any(taxa_sums(rp_lateveg_ps) == 0)
# 4.) Rhizoplane: early reproductive
rp_earrep_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "early_reproductive")
any(taxa_sums(rp_earrep_ps) == 0)
# 5.) Rhizoplane: late reproductive
rp_laterep_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "late_reproductive")
any(taxa_sums(rp_laterep_ps) == 0)



##################################################################################
### Indicator Species Analysis (IndVal.g) - Grouping by Growth Stage Aggregate ###
##################################################################################
library(indicspecies)


### I. RHIZOPLANE 

### Calculate indVal on filtered data ###

# Filter ASVs that had greater than 100 reads (remove Taxa with less than 100 sequences/reads)
rplane_ps_100 <- filter_taxa(rplane_ps, function(x) sum(x) > 100, prune = TRUE)
rplane_ps_100 # 40690 taxa and 133 samples
# Make data frame and transpose to have sample as row names and ASV ID as column names 
asv.rplane.100 <- data.frame(t(otu_table(rplane_ps_100)), check.names = F)
# Subset rhizoplane metadata
rplane_metadata <- data.frame(sample_data(rplane_ps))
# make vector for groups
rplane_stage_aggre <- c(rplane_metadata$growthstage_aggregate_categorical)
# Multi-level pattern analysis: to determine the indicator taxa
set.seed(13)
indic_rplane_stage_aggre_100 <- multipatt(asv.rplane.100, 
                                          rplane_stage_aggre, 
                                          duleg=TRUE, control=how(nperm=999))
# To print only the significant indicator taxa in all growth stages and their Indicator Values:
summary(indic_rplane_stage_aggre_100)
# To extract the complete multipatt results in all growth stages (not only the significant but also the non significant and their Indicator Values):
rp.all_indval.1 <- data.frame(
  asvID = rownames(indic_rplane_stage_aggre_100$sign),
  IndVal.Stat = indic_rplane_stage_aggre_100$sign$stat,
  IndVal.p_value = indic_rplane_stage_aggre_100$sign$p.value,
  Stage = colnames(indic_rplane_stage_aggre_100$comb)[indic_rplane_stage_aggre_100$sign$index])
# Extract multipatt full results in a long format:
# Calculate the Indicator Value score
rp.indval_all <- sqrt(indic_rplane_stage_aggre_100$A * indic_rplane_stage_aggre_100$B)
rp.indval_df <- as.data.frame(rp.indval_all)
rp.indval_df$asvID <- rownames(rp.indval_df)
rp.indval_long <- rp.indval_df %>%
  tidyr::pivot_longer(cols = -asvID, names_to = "Stage", values_to = "IndVal.Stat")
# Add p-value column to the long data frame
rp.indval_long <- rp.indval_long %>%
  left_join(rp.all_indval.1 %>%
              select(asvID, Stage, IndVal.p_value), by = c("asvID", "Stage"))
# Separate by growth stage:
rplane.seedling.indval <- rp.indval_long %>%
  filter(Stage == "seedling")
rplane.early.veget.indval <- rp.indval_long %>%
  filter(Stage == "early_vegetative")
rplane.late.veget.indval <- rp.indval_long %>%
  filter(Stage == "late_vegetative")
rplane.early.reprod.indval <- rp.indval_long %>%
  filter(Stage == "early_reproductive")
rplane.late.reprod.indval <- rp.indval_long %>%
  filter(Stage == "late_reproductive")
# Save the significant only indicator taxa in the computer:
#sink(file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizoplane_indicspec_100.csv")
#summary(indic_rplane_stage_aggre_100) 
#sink()

# Note: Tidy up the data on the Excel, use this code: =TEXTSPLIT(TRIM(A2)," ") in the new column to split between asvID, stage_aggregate, IndVal.Stat, IndVal.p_value, and indVal.signif 
# Read back the edited Indicator ASVs by stage aggregate:
indic_rplane_stage_aggre_100.tdy <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizoplane_indicspec_100_edited.csv")
dim(indic_rplane_stage_aggre_100.tdy) # 905 taxa 
anyDuplicated(indic_rplane_stage_aggre_100.tdy$asvID) == 0 #TRUE: all taxa are unique
# Check how many indicator ASVs in each stage
indic_rplane_stage_aggre_100.tdy$asvID <- as.character(indic_rplane_stage_aggre_100.tdy$asvID)
table(indic_rplane_stage_aggre_100.tdy$stage_aggregate)
# early reproductive   early vegetative  late reproductive    late vegetative           seedling 
#                118                 134                178               293                182


### II. RHIZOSPHERE 

### Calculate indVal on filtered data

# Filter ASVs that had greater than 100 reads
rsphere_ps_100 <- filter_taxa(rsphere_ps, function(x) sum(x) > 100, prune = TRUE)
rsphere_ps_100 # 70986 taxa and 136 samples
# Make data frame and transpose to have sample as row names and ASV ID as column names 
asv.rsphere.100 <- data.frame(t(otu_table(rsphere_ps_100)), check.names = F) 
# Subset rhizoplane metadata
rsphere_metadata <- data.frame(sample_data(rsphere_ps))
# make vector for groups
rsphere_stage_aggre <- c(rsphere_metadata$growthstage_aggregate_categorical)
# Multi-level pattern analysis: to determine the indicator taxa
set.seed(13)
indic_rsphere_stage_aggre_100 <- multipatt(asv.rsphere.100, 
                                           rsphere_stage_aggre, 
                                           duleg=TRUE, control=how(nperm=999))
# To print only the significant indicator taxa in all growth stages and their Indicator Values:
summary(indic_rsphere_stage_aggre_100)
# To extract the complete multipatt results in all growth stages (not only the significant but also the non significant and their Indicator Values):
rs.all_indval.1 <- data.frame(
  asvID = rownames(indic_rsphere_stage_aggre_100$sign),
  IndVal.Stat = indic_rsphere_stage_aggre_100$sign$stat,
  IndVal.p_value = indic_rsphere_stage_aggre_100$sign$p.value,
  Stage = colnames(indic_rsphere_stage_aggre_100$comb)[
    indic_rsphere_stage_aggre_100$sign$index])
# Extract multipatt full results in a long format
# Calculate the Indicator Value score
rs.indval_all <- sqrt(indic_rsphere_stage_aggre_100$A * indic_rsphere_stage_aggre_100$B)
rs.indval_df <- as.data.frame(rs.indval_all)
rs.indval_df$asvID <- rownames(rs.indval_df)
rs.indval_long <- rs.indval_df %>%
  tidyr::pivot_longer(cols = -asvID, names_to = "Stage", values_to = "IndVal.Stat")
# Add p-value column to the long data frame
rs.indval_long <- rs.indval_long %>%
  left_join(rs.all_indval.1 %>%
              select(asvID, Stage, IndVal.p_value), by = c("asvID", "Stage"))
# Separate by growth stage:
rsphere.seedling.indval <- rs.indval_long %>%
  filter(Stage == "seedling")
rsphere.early.veget.indval <- rs.indval_long %>%
  filter(Stage == "early_vegetative")
rsphere.late.veget.indval <- rs.indval_long %>%
  filter(Stage == "late_vegetative")
rsphere.early.reprod.indval <- rs.indval_long %>%
  filter(Stage == "early_reproductive")
rsphere.late.reprod.indval <- rs.indval_long %>%
  filter(Stage == "late_reproductive")
# Save the significant only indicator taxa in the computer:
#sink(file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizosphere_indicspec_100.csv")
#summary(indic_rsphere_stage_aggre_100) 
#sink()

# Note: Tidy up the data on the Excel, use this code: =TEXTSPLIT(TRIM(A2)," ") in the new column to split between asvID, stage_aggregate, IndVal.Stat, IndVal.p_value, and indVal.signif 
# Read back the edited Indicator ASVs by stage aggregate
indic_rsphere_stage_aggre_100.tdy <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizosphere_indicspec_100_edited.csv")
dim(indic_rsphere_stage_aggre_100.tdy) # 872 taxa 
anyDuplicated(indic_rsphere_stage_aggre_100.tdy$asvID) == 0 #TRUE: all taxa are unique
# Check how many indicator ASVs in each stage
indic_rsphere_stage_aggre_100.tdy$asvID <- as.character(indic_rsphere_stage_aggre_100.tdy$asvID)
table(indic_rsphere_stage_aggre_100.tdy$stage_aggregate)
# early reproductive   early vegetative  late reproductive    late vegetative           seedling 
#                100                210                307               86                169




#################################################################################
### Sloan Neutral Community Model (SNCM) - Separate by Growth Stage Aggregate ###
#################################################################################
library(indicspecies)


### I. RHIZOPLANE

# Complete multirarefied rhizoplane phyloseq object
rplane_ps
sample_data(rplane_ps)
# Separate by growth stage aggregate:
rplane_1_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "seedling")
any(taxa_sums(rplane_1_ps) == 0) # 15320 taxa and 30 samples
rplane_2_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "early_vegetative")
any(taxa_sums(rplane_2_ps) == 0) # 14370 taxa and 24 samples
rplane_3_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "late_vegetative")
any(taxa_sums(rplane_3_ps) == 0) # 19387 taxa and 26 samples
rplane_4_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "early_reproductive")
any(taxa_sums(rplane_4_ps) == 0) # 17681 taxa and 29 samples
rplane_5_ps <- rplane_ps %>% ps_filter(growthstage_aggregate_categorical == "late_reproductive")
any(taxa_sums(rplane_5_ps) == 0) # 13081 taxa and 24 samples

# Calculating Mean Relative Abundance and Occupancy; and fit the Sloan Neutral Community Model in each growth stage:

# 1.) Seedling

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the seedling phyloseq object
rplane_1_metadata <- data.frame(sample_data(rplane_1_ps))
sort(colSums(otu_table(rplane_1_ps), na.rm = FALSE, dims = 1), decreasing = T)
# compute relative abundance
otu_table(rplane_1_ps)
ps_1_rel <- transform_sample_counts(
  rplane_1_ps,
  function(x) x / sum(x))
otu_table(ps_1_rel)
# Compute mean relative abundance
rplane_1_asv <- otu_table(rplane_1_ps)
rplane_1_asv_rel <- apply(decostand(rplane_1_asv, method="total", MARGIN=2),1, mean)
rplane_1_asv_rel <- data.frame(rplane_1_asv_rel)
# Compute log10(mean relative abundance)
rplane_1_asv_rel$log_mean <- log10(rplane_1_asv_rel$rplane_1_asv_rel)
rplane_1_asv_rel$asvID <- rownames(rplane_1_asv_rel)
colnames(rplane_1_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rplane_1_asv[rplane_1_asv > 0] <- 1
rplane_1_asv_rel$occupancy <- rowSums(rplane_1_asv)/ncol(rplane_1_asv)
dim(rplane_1_asv_rel) # 15320 
# b) Fit the Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
spp_1 <- data.frame(t(otu_table(rplane_1_ps)), check.names = F)
taxon_1 <- data.frame(tax_table(rplane_1_ps))
# Fit the Sloan Neutral Community Model (SNCM) on seedling
set.seed(13)
rplane_stage1_nm <- fit_sncm(spp_1, pool=NULL, taxon_1)
rplane_stage1_nm$fitstats
# Make the data frame of the SNCM output
rplane_stage1_nm.df <- as.data.frame(rplane_stage1_nm)
arrange(rplane_stage1_nm.df, desc(predictions.freq))
rplane_stage1_nm.df$predictions.fit_class <- as.factor(rplane_stage1_nm.df$predictions.fit_class)
# Indicate taxa that are above CI as TRUE
rplane_stage1_nm.df$above_CI <- rplane_stage1_nm.df$predictions.freq > rplane_stage1_nm.df$predictions.pred.upr
dim(rplane_stage1_nm.df) # 15320 27
# Calculate the percentage of above and below the CI of the model
rplane.stage1.above.pred = sum(rplane_stage1_nm.df$predictions.freq > (rplane_stage1_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rplane_stage1_nm.df)  # fraction of OTUs above prediction # 2.37 %
rplanr.stage1.below.pred = sum(rplane_stage1_nm.df$predictions.freq < (rplane_stage1_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rplane_stage1_nm.df)  # fraction of OTUs below prediction # 1.57 %
100 - (rplane.stage1.above.pred*100) - (rplanr.stage1.below.pred*100)# 96.06%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rplane_stage1_nm.df) == rownames(rplane_1_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rplane_stage1_nm.df), rownames(rplane_1_asv_rel)) # not necessarily the same order
all(rownames(rplane_stage1_nm.df) %in% rownames(rplane_1_asv_rel))  # Should be TRUE
all(rownames(rplane_1_asv_rel) %in% rownames(rplane_stage1_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rplane_stage1_nm.df <- rplane_stage1_nm.df[rownames(rplane_1_asv_rel), , drop = FALSE] # re-order
rplane_stage1_nm.rel <- cbind(rplane_stage1_nm.df, rplane_1_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rplane_stage1_nm.rel[sapply(rplane_stage1_nm.rel, is.character)] <- lapply(rplane_stage1_nm.rel[sapply(rplane_stage1_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Seedling Indicator Taxa Data Frame
# Subset indicator taxa data frame to only seedling indicators
rplane_seedling_indic <- indic_rplane_stage_aggre_100.tdy[indic_rplane_stage_aggre_100.tdy$stage_aggregate == "seedling", ]
dim(rplane_seedling_indic) # 182 seedling indicator taxa
# Check if all seedling indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rplane_stage1_nm.rel <- rownames_to_column(rplane_stage1_nm.rel, var = "asvID")
sum(rplane_seedling_indic$asvID %in% 
      rplane_stage1_nm.rel$asvID) == nrow(rplane_seedling_indic) # TRUE → all seedling indicator ASVs from the "rplane_seedling_indic" are present in the "rplane_stage1_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rplane_stage1_nm.rel.indic <- merge.data.frame(rplane_stage1_nm.rel, 
                                               rplane_seedling_indic, by="asvID", all.x = T)
# Check
sum(rplane_seedling_indic$asvID %in% rplane_stage1_nm.rel$asvID) # n =182, This gives the number of ASVs from "rplane_seedling_indic" that exist in "rplane_stage1_nm.rel".
sum(!rplane_seedling_indic$asvID %in% rplane_stage1_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Seedling
rplane_stage1_nm.rel.indic2 <- rplane_stage1_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rplane.seedling.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rplane_stage1_nm.rel.indic$asvID, rplane_stage1_nm.rel.indic2$asvID)
# Check if all seedling indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rplane_seedling_indic$asvID %in% rplane_stage1_nm.rel.indic2$asvID)
sum(rplane_seedling_indic$asvID %in% 
      rplane_stage1_nm.rel.indic2$asvID) == nrow(rplane_seedling_indic) # TRUE → all seedling indicator ASVs from the "rplane_seedling_indic" are present in the "rplane_stage1_nm.rel"
# Adding stage indicator column (TRUE or FALSE)
rplane_stage1_nm.rel.indic2 <- rplane_stage1_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
dim(rplane_stage1_nm.rel.indic2) #15320 
# Save the complete Seedling taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rplane_stage1_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Seedling_complete_SNCM_Indic.csv")
# Subset taxa that are above CI = deterministic = hypothetically selected by plant
rp.above.stage1 <- subset(rplane_stage1_nm.rel.indic2, above_CI == TRUE)
dim(rp.above.stage1) # 363 taxa are above CI
# e) Identify the overlapping or shared taxa between Indicator Taxa and Above Predicted Taxa in Seedling
# Here we call these taxa as Host-Selected Taxa:
# Subset the data that meet these conditions: above_CI = TRUE ; stage_aggregate = seedling ; stage_indicator = TRUE
rplane_stage1_nm.rel.indic2 <- column_to_rownames(rplane_stage1_nm.rel.indic2, var = "asvID")
rp.indicsp.sncm.overlap.seedling <- subset(rplane_stage1_nm.rel.indic2,
                                           above_CI == TRUE &
                                             stage_aggregate == "seedling" &
                                             stage_indicator == TRUE)
dim(rp.indicsp.sncm.overlap.seedling) # 30 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rp.indicsp.sncm.overlap.seedling, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizoplane_selected_indicator_seedling.csv")

# 2.) Early Vegetative

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the early vegetative phyloseq object
rplane_2_metadata <- data.frame(sample_data(rplane_2_ps))
# Compute mean relative abundance
rplane_2_asv <- otu_table(rplane_2_ps)
rplane_2_asv_rel <- apply(decostand(rplane_2_asv, method="total", MARGIN=2),1, mean)
rplane_2_asv_rel <- data.frame(rplane_2_asv_rel)
# Compute log10(mean relative abundance)
rplane_2_asv_rel$log_mean <- log10(rplane_2_asv_rel$rplane_2_asv_rel)
rplane_2_asv_rel$asvID <- rownames(rplane_2_asv_rel)
colnames(rplane_2_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rplane_2_asv[rplane_2_asv > 0] <- 1
rplane_2_asv_rel$occupancy <- rowSums(rplane_2_asv)/ncol(rplane_2_asv)
dim(rplane_2_asv_rel) # 14370 
# b) Fit the Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
spp_2 <-data.frame(t(otu_table(rplane_2_ps)), check.names = F)
taxon_2 <- data.frame(tax_table(rplane_2_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early vegetative
set.seed(13)
rplane_stage2_nm = fit_sncm(spp_2, pool=NULL, taxon_2)
rplane_stage2_nm$fitstats
# Make the data frame of the SNCM output
rplane_stage2_nm.df <- as.data.frame(rplane_stage2_nm)
arrange(rplane_stage2_nm.df, desc(predictions.freq))
rplane_stage2_nm.df$predictions.fit_class <- as.factor(rplane_stage2_nm.df$predictions.fit_class)
# Indicate taxa that are above CI as TRUE
rplane_stage2_nm.df$above_CI <- rplane_stage2_nm.df$predictions.freq > rplane_stage2_nm.df$predictions.pred.upr
dim(rplane_stage2_nm.df) # 14370 27
# Calculate the percentage of above and below the CI of the model
rplane.stage2.above.pred = sum(rplane_stage2_nm.df$predictions.freq > (rplane_stage2_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rplane_stage2_nm.df)  # fraction of OTUs above prediction # 1.88 %
rplanr.stage2.below.pred = sum(rplane_stage2_nm.df$predictions.freq < (rplane_stage2_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rplane_stage2_nm.df)  # fraction of OTUs below prediction # 1.61 %
100 - (rplane.stage2.above.pred*100) - (rplanr.stage2.below.pred*100)# 96.51%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rplane_stage2_nm.df) == rownames(rplane_2_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rplane_stage2_nm.df), rownames(rplane_2_asv_rel)) # not necessarily the same order
all(rownames(rplane_stage2_nm.df) %in% rownames(rplane_2_asv_rel))  # Should be TRUE
all(rownames(rplane_2_asv_rel) %in% rownames(rplane_stage2_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rplane_stage2_nm.df <- rplane_stage2_nm.df[rownames(rplane_2_asv_rel), , drop = FALSE] # re-order
rplane_stage2_nm.rel <- cbind(rplane_stage2_nm.df, rplane_2_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rplane_stage2_nm.rel[sapply(rplane_stage2_nm.rel, is.character)] <- lapply(rplane_stage2_nm.rel[sapply(rplane_stage2_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Early Vegetative Indicator Taxa Data Frame
# Subset indicator taxa data frame to only early vegetative indicators
rplane_early.veget_indic <- indic_rplane_stage_aggre_100.tdy[indic_rplane_stage_aggre_100.tdy$stage_aggregate == "early_vegetative", ]
dim(rplane_early.veget_indic) # 134 early vegetative indicator taxa
# Check if all early vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rplane_stage2_nm.rel <- rownames_to_column(rplane_stage2_nm.rel, var = "asvID")
sum(rplane_early.veget_indic$asvID %in% rplane_stage2_nm.rel$asvID) == nrow(rplane_early.veget_indic) # TRUE → all early vegetative indicator ASVs from the "rplane_early.veget_indic" are present in the "rplane_stage2_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rplane_stage2_nm.rel.indic <- merge.data.frame(rplane_stage2_nm.rel, rplane_early.veget_indic, by="asvID", all.x = T)
# Check
sum(rplane_early.veget_indic$asvID %in% rplane_stage2_nm.rel$asvID) # n = 134, This gives the number of ASVs from "rplane_early.veget_indic" that exist in "rplane_stage2_nm.rel".
sum(!rplane_early.veget_indic$asvID %in% rplane_stage2_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Early Vegetative
rplane_stage2_nm.rel.indic2 <- rplane_stage2_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rplane.early.veget.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rplane_stage2_nm.rel.indic$asvID, rplane_stage2_nm.rel.indic2$asvID)
# Check if all early vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rplane_early.veget_indic$asvID %in% rplane_stage2_nm.rel.indic2$asvID)
sum(rplane_early.veget_indic$asvID %in% 
      rplane_stage2_nm.rel.indic2$asvID) == nrow(rplane_early.veget_indic) # TRUE 
# Adding stage indicator column (TRUE or FALSE)
rplane_stage2_nm.rel.indic2 <- rplane_stage2_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
dim(rplane_stage2_nm.rel.indic2) # 14370 
# Adding stage indicator column (TRUE or FALSE)
rplane_stage2_nm.rel.indic2 <- rplane_stage2_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Early Vegetative taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rplane_stage2_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Early_Vegetative_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rp.above.stage2 <- subset(rplane_stage2_nm.rel.indic2, above_CI == TRUE)
dim(rp.above.stage2) # 270 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Early vegetative
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = early_vegetative ; stage_indicator = TRUE
# Here we call these taxa as Host-Selected Taxa:
rplane_stage2_nm.rel.indic2 <- column_to_rownames(rplane_stage2_nm.rel.indic2, var = "asvID")
rp.indicsp.sncm.overlap.early.veget <- subset(rplane_stage2_nm.rel.indic2,
                                              above_CI == TRUE &
                                                stage_aggregate == "early_vegetative" &
                                                stage_indicator == TRUE)
dim(rp.indicsp.sncm.overlap.early.veget) # 25 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rp.indicsp.sncm.overlap.early.veget, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizoplane_selected_indicator_early.vegetative.csv")

# 3.) Late Vegetative

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the late vegetative phyloseq object
rplane_3_metadata <- data.frame(sample_data(rplane_3_ps))
# Compute mean relative abundance
rplane_3_asv <- otu_table(rplane_3_ps)
rplane_3_asv_rel <- apply(decostand(rplane_3_asv, method="total", MARGIN=2),1, mean)
rplane_3_asv_rel <- data.frame(rplane_3_asv_rel)
# Compute log10(mean relative abundance)
rplane_3_asv_rel$log_mean <- log10(rplane_3_asv_rel$rplane_3_asv_rel)
rplane_3_asv_rel$asvID <- rownames(rplane_3_asv_rel)
colnames(rplane_3_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rplane_3_asv[rplane_3_asv > 0] <- 1
rplane_3_asv_rel$occupancy <- rowSums(rplane_3_asv)/ncol(rplane_3_asv)
dim(rplane_3_asv_rel) # 19387
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
spp_3 <-data.frame(t(otu_table(rplane_3_ps)), check.names = F)
taxon_3 <- data.frame(tax_table(rplane_3_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early vegetative
set.seed(13)
rplane_stage3_nm = fit_sncm(spp_3, pool=NULL, taxon_3)
rplane_stage3_nm$fitstats
# Make the data frame of the SNCM output
rplane_stage3_nm.df <- as.data.frame(rplane_stage3_nm)
arrange(rplane_stage3_nm.df, desc(predictions.freq))
rplane_stage3_nm.df$predictions.fit_class <- as.factor(rplane_stage3_nm.df$predictions.fit_class)
# Indicate taxa that are above CI as TRUE
rplane_stage3_nm.df$above_CI <- rplane_stage3_nm.df$predictions.freq > rplane_stage3_nm.df$predictions.pred.upr
dim(rplane_stage3_nm.df) # 19387 27
# Calculate the percentage of above and below the CI of the model
rplane.stage3.above.pred = sum(rplane_stage3_nm.df$predictions.freq > (rplane_stage3_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rplane_stage3_nm.df)  # fraction of OTUs above prediction # 1.46 %
rplanr.stage3.below.pred = sum(rplane_stage3_nm.df$predictions.freq < (rplane_stage3_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rplane_stage3_nm.df)  # fraction of OTUs below prediction # 1.34 %
100 - (rplane.stage3.above.pred*100) - (rplanr.stage3.below.pred*100)# 97.19%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rplane_stage3_nm.df) == rownames(rplane_3_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rplane_stage3_nm.df), rownames(rplane_3_asv_rel)) # not necessarily the same order
all(rownames(rplane_stage3_nm.df) %in% rownames(rplane_3_asv_rel))  # Should be TRUE
all(rownames(rplane_3_asv_rel) %in% rownames(rplane_stage3_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rplane_stage3_nm.df <- rplane_stage3_nm.df[rownames(rplane_3_asv_rel), , drop = FALSE] # re-order
rplane_stage3_nm.rel <- cbind(rplane_stage3_nm.df, rplane_3_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rplane_stage3_nm.rel[sapply(rplane_stage3_nm.rel, is.character)] <- lapply(rplane_stage3_nm.rel[sapply(rplane_stage3_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Late Vegetative Indicator Taxa Data Frame
# Subset indicator taxa data frame to only late vegetative indicators
rplane_late.veget_indic <- indic_rplane_stage_aggre_100.tdy[indic_rplane_stage_aggre_100.tdy$stage_aggregate == "late_vegetative", ]
dim(rplane_late.veget_indic) # 293 late vegetative indicator taxa
# Check if all late vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rplane_stage3_nm.rel <- rownames_to_column(rplane_stage3_nm.rel, var = "asvID")
sum(rplane_late.veget_indic$asvID %in% rplane_stage3_nm.rel$asvID) == nrow(rplane_late.veget_indic) # TRUE → all late vegetative indicator ASVs from the "rplane_late.veget_indic" are present in the "rplane_stage3_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rplane_stage3_nm.rel.indic <- merge.data.frame(rplane_stage3_nm.rel, rplane_late.veget_indic, by="asvID", all.x = T)
# Check
sum(rplane_late.veget_indic$asvID %in% rplane_stage3_nm.rel$asvID) # n = 293, This gives the number of ASVs from "rplane_late.veget_indic" that exist in "rplane_stage3_nm.rel".
sum(!rplane_late.veget_indic$asvID %in% rplane_stage3_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Late Vegetative
rplane_stage3_nm.rel.indic2 <- rplane_stage3_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rplane.late.veget.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rplane_stage3_nm.rel.indic$asvID, rplane_stage3_nm.rel.indic2$asvID)
# Check if all late vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rplane_late.veget_indic$asvID %in% rplane_stage3_nm.rel.indic2$asvID)
sum(rplane_late.veget_indic$asvID %in% 
      rplane_stage3_nm.rel.indic2$asvID) == nrow(rplane_late.veget_indic) # TRUE 
# Adding stage indicator column (TRUE or FALSE)
rplane_stage3_nm.rel.indic2 <- rplane_stage3_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
dim(rplane_stage3_nm.rel.indic2) # 19387
# Save the complete Late Vegetative taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rplane_stage3_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Late_Vegetative_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rp.above.stage3 <- subset(rplane_stage3_nm.rel.indic2, above_CI == TRUE)
dim(rp.above.stage3) # 284 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Late vegetative
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = late_vegetative ; stage_indicator = TRUE
# Here we call these taxa as Host-Selected Taxa:
rplane_stage3_nm.rel.indic2 <- column_to_rownames(rplane_stage3_nm.rel.indic2, var = "asvID")
rp.indicsp.sncm.overlap.late.veget <- subset(rplane_stage3_nm.rel.indic2,
                                             above_CI == TRUE &
                                               stage_aggregate == "late_vegetative" &
                                               stage_indicator == TRUE)
dim(rp.indicsp.sncm.overlap.late.veget) # 54 Host-Selected Taxa:
# Save the Host-Selected Taxa:
#write.csv(rp.indicsp.sncm.overlap.late.veget, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizoplane_selected_indicator_late.vegetative.csv")

# 4.) Early Reproductive

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the early reproductive phyloseq object
rplane_4_metadata <- data.frame(sample_data(rplane_4_ps))
# Compute mean relative abundance
rplane_4_asv <- otu_table(rplane_4_ps)
rplane_4_asv_rel <- apply(decostand(rplane_4_asv, method="total", MARGIN=2),1, mean)
rplane_4_asv_rel <- data.frame(rplane_4_asv_rel)
# Compute log10(mean relative abundance)
rplane_4_asv_rel$log_mean <- log10(rplane_4_asv_rel$rplane_4_asv_rel)
rplane_4_asv_rel$asvID <- rownames(rplane_4_asv_rel)
colnames(rplane_4_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rplane_4_asv[rplane_4_asv > 0] <- 1
rplane_4_asv_rel$occupancy <- rowSums(rplane_4_asv)/ncol(rplane_4_asv)
dim(rplane_4_asv_rel) # 17681 
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
spp_4 <-data.frame(t(otu_table(rplane_4_ps)), check.names = F)
taxon_4 <- data.frame(tax_table(rplane_4_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early reproductive
set.seed(13)
rplane_stage4_nm = fit_sncm(spp_4, pool=NULL, taxon_4)
rplane_stage4_nm$fitstats
# Make the data frame of the SNCM output
rplane_stage4_nm.df <- as.data.frame(rplane_stage4_nm)
arrange(rplane_stage4_nm.df, desc(predictions.freq))
rplane_stage4_nm.df$predictions.fit_class <- as.factor(rplane_stage4_nm.df$predictions.fit_class)
# Indicate taxa that are above CI as TRUE
rplane_stage4_nm.df$above_CI <- rplane_stage4_nm.df$predictions.freq > rplane_stage4_nm.df$predictions.pred.upr
dim(rplane_stage4_nm.df) # 17681 27
# Calculate the percentage of above and below the CI of the model
rplane.stage4.above.pred = sum(rplane_stage4_nm.df$predictions.freq > (rplane_stage4_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rplane_stage4_nm.df)  # fraction of OTUs above prediction # 1.85 %
rplanr.stage4.below.pred = sum(rplane_stage4_nm.df$predictions.freq < (rplane_stage4_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rplane_stage4_nm.df)  # fraction of OTUs below prediction # 1.33 %
100 - (rplane.stage4.above.pred*100) - (rplanr.stage4.below.pred*100)# 96.82%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rplane_stage4_nm.df) == rownames(rplane_4_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rplane_stage4_nm.df), rownames(rplane_4_asv_rel)) # not necessarily the same order
all(rownames(rplane_stage4_nm.df) %in% rownames(rplane_4_asv_rel))  # Should be TRUE
all(rownames(rplane_4_asv_rel) %in% rownames(rplane_stage4_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rplane_stage4_nm.df <- rplane_stage4_nm.df[rownames(rplane_4_asv_rel), , drop = FALSE] # re-order
rplane_stage4_nm.rel <- cbind(rplane_stage4_nm.df, rplane_4_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rplane_stage4_nm.rel[sapply(rplane_stage4_nm.rel, is.character)] <- lapply(rplane_stage4_nm.rel[sapply(rplane_stage4_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Early Reproductive Indicator Taxa Data Frame
# Subset indicator taxa data frame to only early reproductive indicators
rplane_early.reprod_indic <- indic_rplane_stage_aggre_100.tdy[indic_rplane_stage_aggre_100.tdy$stage_aggregate == "early_reproductive", ]
dim(rplane_early.reprod_indic) # 118 early reproductive indicator taxa
# Check if all early reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rplane_stage4_nm.rel <- rownames_to_column(rplane_stage4_nm.rel, var = "asvID")
sum(rplane_early.reprod_indic$asvID %in% rplane_stage4_nm.rel$asvID) == nrow(rplane_early.reprod_indic) # TRUE → all early reproductive indicator ASVs from the "rplane_early.reprod_indic" are present in the "rplane_stage4_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rplane_stage4_nm.rel.indic <- merge.data.frame(rplane_stage4_nm.rel, rplane_early.reprod_indic, by="asvID", all.x = T)
# Check
sum(rplane_early.reprod_indic$asvID %in% rplane_stage4_nm.rel$asvID) # n = 118, This gives the number of ASVs from "rplane_early.reprod_indic" that exist in "rplane_stage4_nm.rel".
sum(!rplane_early.reprod_indic$asvID %in% rplane_stage4_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Early Reproductive
rplane_stage4_nm.rel.indic2 <- rplane_stage4_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rplane.early.reprod.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rplane_stage4_nm.rel.indic$asvID, rplane_stage4_nm.rel.indic2$asvID)
# Check if all early reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rplane_early.reprod_indic$asvID %in% rplane_stage4_nm.rel.indic2$asvID)
sum(rplane_early.reprod_indic$asvID %in% 
      rplane_stage4_nm.rel.indic2$asvID) == nrow(rplane_early.reprod_indic) # TRUE 
# Adding stage indicator column (TRUE or FALSE)
rplane_stage4_nm.rel.indic2 <- rplane_stage4_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
dim(rplane_stage4_nm.rel.indic2) # 17681
# Save the complete Early Reproductive taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rplane_stage4_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Early_Reproductive_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rp.above.stage4 <- subset(rplane_stage4_nm.rel.indic2, above_CI == TRUE)
dim(rp.above.stage4) # 327 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Early reproductive
# Subset the data that meet these conditions: above_CI = TRUE ; stage_aggregate = early_reproductive ; stage_indicator = TRUE
# Here we call these taxa as Host-Selected Taxa:
rplane_stage4_nm.rel.indic2 <- column_to_rownames(rplane_stage4_nm.rel.indic2, var = "asvID")
rp.indicsp.sncm.overlap.early.reprod <- subset(rplane_stage4_nm.rel.indic2,
                                               above_CI == TRUE &
                                                 stage_aggregate == "early_reproductive" &
                                                 stage_indicator == TRUE)
dim(rp.indicsp.sncm.overlap.early.reprod) # 26 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rp.indicsp.sncm.overlap.early.reprod, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizoplane_selected_indicator_early.reproductive.csv")

# 5.) Late Reproductive

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the late reproductive phyloseq object
rplane_5_metadata <- data.frame(sample_data(rplane_5_ps))
# Compute mean relative abundance
rplane_5_asv <- otu_table(rplane_5_ps)
rplane_5_asv_rel <- apply(decostand(rplane_5_asv, method="total", MARGIN=2),1, mean)
rplane_5_asv_rel <- data.frame(rplane_5_asv_rel)
# Compute log10(mean relative abundance)
rplane_5_asv_rel$log_mean <- log10(rplane_5_asv_rel$rplane_5_asv_rel)
rplane_5_asv_rel$asvID <- rownames(rplane_5_asv_rel)
colnames(rplane_5_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rplane_5_asv[rplane_5_asv > 0] <- 1
rplane_5_asv_rel$occupancy <- rowSums(rplane_5_asv)/ncol(rplane_5_asv)
dim(rplane_5_asv_rel) # 13081
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
spp_5 <-data.frame(t(otu_table(rplane_5_ps)), check.names = F)
taxon_5 <- data.frame(tax_table(rplane_5_ps))
# Fit the Sloan Neutral Community Model (SNCM) on late reproductive
set.seed(13)
rplane_stage5_nm <- fit_sncm(spp_5, pool=NULL, taxon_5)
rplane_stage5_nm$fitstats
# Make the data frame of the SNCM output
rplane_stage5_nm.df <- as.data.frame(rplane_stage5_nm)
arrange(rplane_stage5_nm.df, desc(predictions.freq))
rplane_stage5_nm.df$predictions.fit_class <- as.factor(rplane_stage5_nm.df$predictions.fit_class)
# Indicate taxa that are above CI as TRUE
rplane_stage5_nm.df$above_CI <- rplane_stage5_nm.df$predictions.freq > rplane_stage5_nm.df$predictions.pred.upr
dim(rplane_stage5_nm.df) # 13081 27
# Calculate the percentage of above and below the CI of the model
rplane.stage5.above.pred = sum(rplane_stage5_nm.df$predictions.freq > (rplane_stage5_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rplane_stage5_nm.df)  # fraction of OTUs above prediction # 1.76 %
rplanr.stage5.below.pred = sum(rplane_stage5_nm.df$predictions.freq < (rplane_stage5_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rplane_stage5_nm.df)  # fraction of OTUs below prediction # 1.39 %
100 - (rplane.stage5.above.pred*100) - (rplanr.stage5.below.pred*100)# 96.84%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rplane_stage5_nm.df) == rownames(rplane_5_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rplane_stage5_nm.df), rownames(rplane_5_asv_rel)) # not necessarily the same order
all(rownames(rplane_stage5_nm.df) %in% rownames(rplane_5_asv_rel))  # Should be TRUE
all(rownames(rplane_5_asv_rel) %in% rownames(rplane_stage5_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rplane_stage5_nm.df <- rplane_stage5_nm.df[rownames(rplane_5_asv_rel), , drop = FALSE] # re-order
rplane_stage5_nm.rel <- cbind(rplane_stage5_nm.df, rplane_5_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rplane_stage5_nm.rel[sapply(rplane_stage5_nm.rel, is.character)] <- lapply(rplane_stage5_nm.rel[sapply(rplane_stage5_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Late Reproductive Indicator Taxa Data Frame
# Subset indicator taxa data frame to only late reproductive indicators
rplane_late.reprod_indic <- indic_rplane_stage_aggre_100.tdy[indic_rplane_stage_aggre_100.tdy$stage_aggregate == "late_reproductive", ]
dim(rplane_late.reprod_indic) # 178 late reproductive indicator taxa
# Check if all late reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rplane_stage5_nm.rel <- rownames_to_column(rplane_stage5_nm.rel, var = "asvID")
sum(rplane_late.reprod_indic$asvID %in% rplane_stage5_nm.rel$asvID) == nrow(rplane_late.reprod_indic) # TRUE → all late reproductive indicator ASVs from the "rplane_late.reprod_indic" are present in the "rplane_stage5_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rplane_stage5_nm.rel.indic <- merge.data.frame(rplane_stage5_nm.rel, rplane_late.reprod_indic, by="asvID", all.x = T)
# Check
sum(rplane_late.reprod_indic$asvID %in% rplane_stage5_nm.rel$asvID) # n = 178, This gives the number of ASVs from "rplane_late.reprod_indic" that exist in "rplane_stage5_nm.rel".
sum(!rplane_late.reprod_indic$asvID %in% rplane_stage5_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Late Reproductive
rplane_stage5_nm.rel.indic2 <- rplane_stage5_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rplane.late.reprod.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rplane_stage5_nm.rel.indic$asvID, rplane_stage5_nm.rel.indic2$asvID)
# Check if all late reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rplane_late.reprod_indic$asvID %in% rplane_stage5_nm.rel.indic2$asvID)
sum(rplane_late.reprod_indic$asvID %in% 
      rplane_stage5_nm.rel.indic2$asvID) == nrow(rplane_late.reprod_indic) # TRUE 
# Adding stage indicator column (TRUE or FALSE)
rplane_stage5_nm.rel.indic2 <- rplane_stage5_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Late Reproductive taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rplane_stage5_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizoplane_Late_Reproductive_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rp.above.stage5 <- subset(rplane_stage5_nm.rel.indic2, above_CI == TRUE)
dim(rp.above.stage5) # 230 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Late reproductive
# Subset the data that meet these conditions: above_CI = TRUE ; stage_aggregate = late_reproductive ; stage_indicator = TRUE
# Here we call these taxa as Host-Selected Taxa:
rplane_stage5_nm.rel.indic2 <- column_to_rownames(rplane_stage5_nm.rel.indic2, var = "asvID")
rp.indicsp.sncm.overlap.late.reprod <- subset(rplane_stage5_nm.rel.indic2,
                                              above_CI == TRUE &
                                                stage_aggregate == "late_reproductive" &
                                                stage_indicator == TRUE)
dim(rp.indicsp.sncm.overlap.late.reprod) # 24 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rp.indicsp.sncm.overlap.late.reprod, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizoplane_selected_indicator_late.reproductive.csv")


### 2. RHIZOSPHERE

# Complete multirarefied rhizoplane phyloseq object
rsphere_ps # 85345 taxa and 136 samples 
sample_data(rsphere_ps)
# Separate by growth stage aggregate
rsphere_1_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "seedling")
any(taxa_sums(rsphere_1_ps) == 0) # 21962 taxa and 30 samples
rsphere_2_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "early_vegetative")
any(taxa_sums(rsphere_2_ps) == 0) # 17980 taxa and 22 samples
rsphere_3_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "late_vegetative")
any(taxa_sums(rsphere_3_ps) == 0) # 23084 taxa and 30 samples
rsphere_4_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "early_reproductive")
any(taxa_sums(rsphere_4_ps) == 0) # 21987 taxa and 29 samples
rsphere_5_ps <- rsphere_ps %>% ps_filter(growthstage_aggregate_categorical == "late_reproductive")
any(taxa_sums(rsphere_5_ps) == 0) # 19549 taxa and 25 samples

# Calculating Mean Relative Abundance and Occupancy; and fit the Sloan Neutral Community Model in each growth stage:

### 1.) Seedling

# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the seedling phyloseq object
rsphere_1_metadata <- data.frame(sample_data(rsphere_1_ps))
# Compute mean relative abundance
rsphere_1_asv <- otu_table(rsphere_1_ps)
rsphere_1_asv_rel <- apply(decostand(rsphere_1_asv, method="total", MARGIN=2),1, mean)
rsphere_1_asv_rel <- data.frame(rsphere_1_asv_rel)
# Compute log10(mean relative abundance)
rsphere_1_asv_rel$log_mean <- log10(rsphere_1_asv_rel$rsphere_1_asv_rel)
rsphere_1_asv_rel$asvID <- rownames(rsphere_1_asv_rel)
colnames(rsphere_1_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rsphere_1_asv[rsphere_1_asv > 0] <- 1
rsphere_1_asv_rel$occupancy <- rowSums(rsphere_1_asv)/ncol(rsphere_1_asv)
dim(rsphere_1_asv_rel) # 21962 
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
rsphere.spp_1 <- data.frame(t(otu_table(rsphere_1_ps)), check.names = F)
rsphere.taxon_1 <- data.frame(tax_table(rsphere_1_ps))
# Fit the Sloan Neutral Community Model (SNCM) on seedling
set.seed(13)
rsphere_stage1_nm <- fit_sncm(rsphere.spp_1, pool=NULL, rsphere.taxon_1)
rsphere_stage1_nm$fitstats
# Make the data frame of the SNCM output
rsphere_stage1_nm.df <- as.data.frame(rsphere_stage1_nm)
arrange(rsphere_stage1_nm.df, desc(predictions.freq))
rsphere_stage1_nm.df$predictions.fit_class <- as.factor(rsphere_stage1_nm.df$predictions.fit_class)
rsphere_stage1_nm.df$predictions.fit_class 
# Indicate taxa that are above CI as TRUE
rsphere_stage1_nm.df$above_CI <- rsphere_stage1_nm.df$predictions.freq > rsphere_stage1_nm.df$predictions.pred.upr
dim(rsphere_stage1_nm.df) # 21962 27
# Calculate the percentage of above and below the CI of the model
rsphere.stage1.above.pred = sum(rsphere_stage1_nm.df$predictions.freq > (rsphere_stage1_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rsphere_stage1_nm.df)  # fraction of OTUs above prediction # 1.55 %
rsphere.stage1.below.pred = sum(rsphere_stage1_nm.df$predictions.freq < (rsphere_stage1_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rsphere_stage1_nm.df)  # fraction of OTUs below prediction # 1.38 %
100 - (rsphere.stage1.above.pred*100) - (rsphere.stage1.below.pred*100)# 97.08%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rsphere_stage1_nm.df) == rownames(rsphere_1_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rsphere_stage1_nm.df), rownames(rsphere_1_asv_rel)) # not necessarily the same order
all(rownames(rsphere_stage1_nm.df) %in% rownames(rsphere_1_asv_rel))  # Should be TRUE
all(rownames(rsphere_1_asv_rel) %in% rownames(rsphere_stage1_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rsphere_stage1_nm.df <- rsphere_stage1_nm.df[rownames(rsphere_1_asv_rel), , drop = FALSE] # re-order
rsphere_stage1_nm.rel <- cbind(rsphere_stage1_nm.df, rsphere_1_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rsphere_stage1_nm.rel[sapply(rsphere_stage1_nm.rel, is.character)] <- lapply(rsphere_stage1_nm.rel[sapply(rsphere_stage1_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Seedling Indicator Taxa Data Frame
# Subset indicator taxa data frame to only seedling indicators
rsphere_seedling_indic <- indic_rsphere_stage_aggre_100.tdy[indic_rsphere_stage_aggre_100.tdy$stage_aggregate == "seedling", ]
dim(rsphere_seedling_indic) # 169 seedling indicator taxa
# Check if all seedling indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rsphere_stage1_nm.rel <- rownames_to_column(rsphere_stage1_nm.rel, var = "asvID")
sum(rsphere_seedling_indic$asvID %in% rsphere_stage1_nm.rel$asvID) == nrow(rsphere_seedling_indic) # TRUE → all seedling indicator ASVs from the "rsphere_seedling_indic" are present in the "rsphere_stage1_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rsphere_stage1_nm.rel.indic <- merge.data.frame(rsphere_stage1_nm.rel, rsphere_seedling_indic, by="asvID", all.x = T)
# Check
sum(rsphere_seedling_indic$asvID %in% rsphere_stage1_nm.rel$asvID) # n = 169, This gives the number of ASVs from "rsphere_seedling_indic" that exist in "rsphere_stage1_nm.rel".
sum(!rsphere_seedling_indic$asvID %in% rsphere_stage1_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Seedling
rsphere_stage1_nm.rel.indic2 <- rsphere_stage1_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rsphere.seedling.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rsphere_stage1_nm.rel.indic$asvID, rsphere_stage1_nm.rel.indic2$asvID)
# Check if all seedling indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rsphere_seedling_indic$asvID %in% rsphere_stage1_nm.rel.indic2$asvID)
sum(rsphere_seedling_indic$asvID %in% 
      rsphere_stage1_nm.rel.indic2$asvID) == nrow(rsphere_seedling_indic) # TRUE 
# Adding stage indicator column (TRUE or FALSE)
rsphere_stage1_nm.rel.indic2 <- rsphere_stage1_nm.rel.indic2 %>%
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Seedling taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rsphere_stage1_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Seedling_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rsphere.above.stage1 <- subset(rsphere_stage1_nm.rel.indic2, above_CI == TRUE)
dim(rsphere.above.stage1) # 340 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Seedling
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = seedling ; stage_indicator = TRUE
rsphere_stage1_nm.rel.indic2 <- column_to_rownames(rsphere_stage1_nm.rel.indic2, var = "asvID")
rsphere.indicsp.sncm.overlap.seedling <- subset(rsphere_stage1_nm.rel.indic2,
                                                above_CI == TRUE &
                                                  stage_aggregate == "seedling" &
                                                  stage_indicator == TRUE)
dim(rsphere.indicsp.sncm.overlap.seedling) # 29 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rsphere.indicsp.sncm.overlap.seedling, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizosphere_selected_indicator_seedling.csv")

### 2.) Early Vegetative
# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the early vegetative phyloseq object
rsphere_2_metadata <- data.frame(sample_data(rsphere_2_ps))
# Compute mean relative abundance
rsphere_2_asv <- otu_table(rsphere_2_ps)
rsphere_2_asv_rel <- apply(decostand(rsphere_2_asv, method="total", MARGIN=2),1, mean)
rsphere_2_asv_rel <- data.frame(rsphere_2_asv_rel)
# Compute log10(mean relative abundance)
rsphere_2_asv_rel$log_mean <- log10(rsphere_2_asv_rel$rsphere_2_asv_rel)
rsphere_2_asv_rel$asvID <- rownames(rsphere_2_asv_rel)
colnames(rsphere_2_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rsphere_2_asv[rsphere_2_asv > 0] <- 1
rsphere_2_asv_rel$occupancy <- rowSums(rsphere_2_asv)/ncol(rsphere_2_asv)
dim(rsphere_2_asv_rel) # 17980 
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
rsphere.spp_2 <-data.frame(t(otu_table(rsphere_2_ps)), check.names = F)
rsphere.taxon_2 <- data.frame(tax_table(rsphere_2_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early vegetative
set.seed(13)
rsphere_stage2_nm = fit_sncm(rsphere.spp_2, pool=NULL, rsphere.taxon_2)
rsphere_stage2_nm$fitstats
# Make the data frame of the SNCM output
rsphere_stage2_nm.df <- as.data.frame(rsphere_stage2_nm)
arrange(rsphere_stage2_nm.df, desc(predictions.freq))
rsphere_stage2_nm.df$predictions.fit_class <- as.factor(rsphere_stage2_nm.df$predictions.fit_class)
rsphere_stage2_nm.df$predictions.fit_class 
# Indicate taxa that are above CI as TRUE
rsphere_stage2_nm.df$above_CI <- rsphere_stage2_nm.df$predictions.freq > rsphere_stage2_nm.df$predictions.pred.upr
dim(rsphere_stage2_nm.df) # 17980 27
# Calculate the percentage of above and below the CI of the model
rsphere.stage2.above.pred = sum(rsphere_stage2_nm.df$predictions.freq > (rsphere_stage2_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rsphere_stage2_nm.df)  # fraction of OTUs above prediction # 1.51 %
rsphere.stage2.below.pred = sum(rsphere_stage2_nm.df$predictions.freq < (rsphere_stage2_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rsphere_stage2_nm.df)  # fraction of OTUs below prediction # 1.39 %
100 - (rsphere.stage2.above.pred*100) - (rsphere.stage2.below.pred*100)# 97.10%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rsphere_stage2_nm.df) == rownames(rsphere_2_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rsphere_stage2_nm.df), rownames(rsphere_2_asv_rel)) # not necessarily the same order
all(rownames(rsphere_stage2_nm.df) %in% rownames(rsphere_2_asv_rel))  # Should be TRUE
all(rownames(rsphere_2_asv_rel) %in% rownames(rsphere_stage2_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rsphere_stage2_nm.df <- rsphere_stage2_nm.df[rownames(rsphere_2_asv_rel), , drop = FALSE] # re-order
rsphere_stage2_nm.rel <- cbind(rsphere_stage2_nm.df, rsphere_2_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rsphere_stage2_nm.rel[sapply(rsphere_stage2_nm.rel, is.character)] <- lapply(rsphere_stage2_nm.rel[sapply(rsphere_stage2_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Early Vegetative Indicator Taxa Data Frame
# Subset indicator taxa data frame to only early vegetative indicators
rsphere_early.veget_indic <- indic_rsphere_stage_aggre_100.tdy[indic_rsphere_stage_aggre_100.tdy$stage_aggregate == "early_vegetative", ]
dim(rsphere_early.veget_indic) # 210 early vegetative indicator taxa
# Check if all early vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rsphere_stage2_nm.rel <- rownames_to_column(rsphere_stage2_nm.rel, var = "asvID")
sum(rsphere_early.veget_indic$asvID %in% rsphere_stage2_nm.rel$asvID) == nrow(rsphere_early.veget_indic) # TRUE → all early vegetative indicator ASVs from the "rsphere_early.veget_indic" are present in the "rsphere_stage2_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rsphere_stage2_nm.rel.indic <- merge.data.frame(rsphere_stage2_nm.rel, rsphere_early.veget_indic, by="asvID", all.x = T)
# Check
sum(rsphere_early.veget_indic$asvID %in% rsphere_stage2_nm.rel$asvID) # n = 210, This gives the number of ASVs from "rsphere_early.veget_indic" that exist in "rsphere_stage2_nm.rel".
sum(!rsphere_early.veget_indic$asvID %in% rsphere_stage2_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Early Vegetative
rsphere_stage2_nm.rel.indic2 <- rsphere_stage2_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rsphere.early.veget.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rsphere_stage2_nm.rel.indic$asvID, rsphere_stage2_nm.rel.indic2$asvID)
# Check if all early vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rsphere_early.veget_indic$asvID %in% rsphere_stage2_nm.rel.indic2$asvID)
sum(rsphere_early.veget_indic$asvID %in% 
      rsphere_stage2_nm.rel.indic2$asvID) == nrow(rsphere_early.veget_indic) #TRUE
# Adding stage indicator column (TRUE or FALSE)
rsphere_stage2_nm.rel.indic2 <- rsphere_stage2_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Early Vegetative taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rsphere_stage2_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Early_Vegetative_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rsphere.above.stage2 <- subset(rsphere_stage2_nm.rel.indic2, above_CI == TRUE)
dim(rsphere.above.stage2) # 271 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Early vegetative
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = early_vegetative ; stage_indicator = TRUE
rsphere_stage2_nm.rel.indic2 <- column_to_rownames(rsphere_stage2_nm.rel.indic2, var = "asvID")
rsphere.indicsp.sncm.overlap.early.veget <- subset(rsphere_stage2_nm.rel.indic2,
                                                   above_CI == TRUE &
                                                     stage_aggregate == "early_vegetative" &
                                                     stage_indicator == TRUE)
rownames(rsphere.indicsp.sncm.overlap.early.veget) # 32 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rsphere.indicsp.sncm.overlap.early.veget, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizosphere_selected_indicator_early.vegetative.csv")

### 3.) Late Vegetative
# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the late vegetative phyloseq object
rsphere_3_metadata <- data.frame(sample_data(rsphere_3_ps))
# Compute mean relative abundance
rsphere_3_asv <- otu_table(rsphere_3_ps)
rsphere_3_asv_rel <- apply(decostand(rsphere_3_asv, method="total", MARGIN=2),1, mean)
rsphere_3_asv_rel <- data.frame(rsphere_3_asv_rel)
# Compute log10(mean relative abundance)
rsphere_3_asv_rel$log_mean <- log10(rsphere_3_asv_rel$rsphere_3_asv_rel)
rsphere_3_asv_rel$asvID <- rownames(rsphere_3_asv_rel)
colnames(rsphere_3_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rsphere_3_asv[rsphere_3_asv > 0] <- 1
rsphere_3_asv_rel$occupancy <- rowSums(rsphere_3_asv)/ncol(rsphere_3_asv)
dim(rsphere_3_asv_rel) # 23084
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
rsphere.spp_3 <-data.frame(t(otu_table(rsphere_3_ps)), check.names = F)
rsphere.taxon_3 <- data.frame(tax_table(rsphere_3_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early vegetative
set.seed(13)
rsphere_stage3_nm = fit_sncm(rsphere.spp_3, pool=NULL, rsphere.taxon_3)
rsphere_stage3_nm$fitstats
# Make the data frame of the SNCM output
rsphere_stage3_nm.df <- as.data.frame(rsphere_stage3_nm)
arrange(rsphere_stage3_nm.df, desc(predictions.freq))
rsphere_stage3_nm.df$predictions.fit_class <- as.factor(rsphere_stage3_nm.df$predictions.fit_class)
rsphere_stage3_nm.df$predictions.fit_class 
# Indicate taxa that are above CI as TRUE
rsphere_stage3_nm.df$above_CI <- rsphere_stage3_nm.df$predictions.freq > rsphere_stage3_nm.df$predictions.pred.upr
dim(rsphere_stage3_nm.df) # 23084 27
# Calculate the percentage of above and below the CI of the model
rsphere.stage3.above.pred = sum(rsphere_stage3_nm.df$predictions.freq > (rsphere_stage3_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rsphere_stage3_nm.df)  # fraction of OTUs above prediction # 1.68 %
rsphere.stage3.below.pred = sum(rsphere_stage3_nm.df$predictions.freq < (rsphere_stage3_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rsphere_stage3_nm.df)  # fraction of OTUs below prediction # 1.37 %
100 - (rsphere.stage3.above.pred*100) - (rsphere.stage3.below.pred*100)# 96.95%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rsphere_stage3_nm.df) == rownames(rsphere_3_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rsphere_stage3_nm.df), rownames(rsphere_3_asv_rel)) # not necessarily the same order
all(rownames(rsphere_stage3_nm.df) %in% rownames(rsphere_3_asv_rel))  # Should be TRUE
all(rownames(rsphere_3_asv_rel) %in% rownames(rsphere_stage3_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rsphere_stage3_nm.df <- rsphere_stage3_nm.df[rownames(rsphere_3_asv_rel), , drop = FALSE] # re-order
rsphere_stage3_nm.rel <- cbind(rsphere_stage3_nm.df, rsphere_3_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rsphere_stage3_nm.rel[sapply(rsphere_stage3_nm.rel, is.character)] <- lapply(rsphere_stage3_nm.rel[sapply(rsphere_stage3_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Late Vegetative Indicator Taxa Data Frame
# Subset indicator taxa data frame to only late vegetative indicators
rsphere_late.veget_indic <- indic_rsphere_stage_aggre_100.tdy[indic_rsphere_stage_aggre_100.tdy$stage_aggregate == "late_vegetative", ]
dim(rsphere_late.veget_indic) # 86 late vegetative indicator taxa
# Check if all late vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rsphere_stage3_nm.rel <- rownames_to_column(rsphere_stage3_nm.rel, var = "asvID")
sum(rsphere_late.veget_indic$asvID %in% rsphere_stage3_nm.rel$asvID) == nrow(rsphere_late.veget_indic) # TRUE → all late vegetative indicator ASVs from the "rsphere_late.veget_indic" are present in the "rsphere_stage3_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rsphere_stage3_nm.rel.indic <- merge.data.frame(rsphere_stage3_nm.rel, rsphere_late.veget_indic, by="asvID", all.x = T)
# Check
sum(rsphere_late.veget_indic$asvID %in% rsphere_stage3_nm.rel$asvID) # n = 86, This gives the number of ASVs from "rsphere_late.veget_indic" that exist in "rsphere_stage3_nm.rel".
sum(!rsphere_late.veget_indic$asvID %in% rsphere_stage3_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Late Vegetative
rsphere_stage3_nm.rel.indic2 <- rsphere_stage3_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rsphere.late.veget.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rsphere_stage3_nm.rel.indic$asvID, rsphere_stage3_nm.rel.indic2$asvID)
# Check if all early vegetative indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rsphere_late.veget_indic$asvID %in% rsphere_stage3_nm.rel.indic2$asvID)
sum(rsphere_late.veget_indic$asvID %in% 
      rsphere_stage3_nm.rel.indic2$asvID) == nrow(rsphere_late.veget_indic) #TRUE
# Adding stage indicator column (TRUE or FALSE)
rsphere_stage3_nm.rel.indic2 <- rsphere_stage3_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Late Vegetative taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rsphere_stage3_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Late_Vegetative_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rsphere.above.stage3 <- subset(rsphere_stage3_nm.rel.indic2, above_CI == TRUE)
dim(rsphere.above.stage3) # 388 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Late vegetative
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = late_vegetative ; stage_indicator = TRUE
rsphere_stage3_nm.rel.indic2 <- column_to_rownames(rsphere_stage3_nm.rel.indic2, var = "asvID")
rsphere.indicsp.sncm.overlap.late.veget <- subset(rsphere_stage3_nm.rel.indic2,
                                                  above_CI == TRUE &
                                                    stage_aggregate == "late_vegetative" &
                                                    stage_indicator == TRUE)
dim(rsphere.indicsp.sncm.overlap.late.veget) # 18 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rsphere.indicsp.sncm.overlap.late.veget, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizosphere_selected_indicator_late.vegetative.csv")

### 4.) Early Reproductive
# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the early reproductive phyloseq object
rsphere_4_metadata <- data.frame(sample_data(rsphere_4_ps))
# Compute mean relative abundance
rsphere_4_asv <- otu_table(rsphere_4_ps)
rsphere_4_asv_rel <- apply(decostand(rsphere_4_asv, method="total", MARGIN=2),1, mean)
rsphere_4_asv_rel <- data.frame(rsphere_4_asv_rel)
# Compute log10(mean relative abundance)
rsphere_4_asv_rel$log_mean <- log10(rsphere_4_asv_rel$rsphere_4_asv_rel)
rsphere_4_asv_rel$asvID <- rownames(rsphere_4_asv_rel)
colnames(rsphere_4_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rsphere_4_asv[rsphere_4_asv > 0] <- 1
rsphere_4_asv_rel$occupancy <- rowSums(rsphere_4_asv)/ncol(rsphere_4_asv)
dim(rsphere_4_asv_rel) # 21987 
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
rsphere.spp_4 <-data.frame(t(otu_table(rsphere_4_ps)), check.names = F)
rsphere.taxon_4 <- data.frame(tax_table(rsphere_4_ps))
# Fit the Sloan Neutral Community Model (SNCM) on early reproductive
set.seed(13)
rsphere_stage4_nm = fit_sncm(rsphere.spp_4, pool=NULL, rsphere.taxon_4)
rsphere_stage4_nm$fitstats
# Make the data frame of the SNCM output
rsphere_stage4_nm.df <- as.data.frame(rsphere_stage4_nm)
arrange(rsphere_stage4_nm.df, desc(predictions.freq))
rsphere_stage4_nm.df$predictions.fit_class <- as.factor(rsphere_stage4_nm.df$predictions.fit_class)
rsphere_stage4_nm.df$predictions.fit_class 
# Indicate taxa that are above CI as TRUE
rsphere_stage4_nm.df$above_CI <- rsphere_stage4_nm.df$predictions.freq > rsphere_stage4_nm.df$predictions.pred.upr
dim(rsphere_stage4_nm.df) # 21987 27
# Calculate the percentage of above and below the CI of the model
rsphere.stage4.above.pred = sum(rsphere_stage4_nm.df$predictions.freq > (rsphere_stage4_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rsphere_stage4_nm.df)  # fraction of OTUs above prediction # 1.56 %
rsphere.stage4.below.pred = sum(rsphere_stage4_nm.df$predictions.freq < (rsphere_stage4_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rsphere_stage4_nm.df)  # fraction of OTUs below prediction # 1.26 %
100 - (rsphere.stage4.above.pred*100) - (rsphere.stage4.below.pred*100)# 97.19%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rsphere_stage4_nm.df) == rownames(rsphere_4_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rsphere_stage4_nm.df), rownames(rsphere_4_asv_rel)) # not necessarily the same order
all(rownames(rsphere_stage4_nm.df) %in% rownames(rsphere_4_asv_rel))  # Should be TRUE
all(rownames(rsphere_4_asv_rel) %in% rownames(rsphere_stage4_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rsphere_stage4_nm.df <- rsphere_stage4_nm.df[rownames(rsphere_4_asv_rel), , drop = FALSE] # re-order
rsphere_stage4_nm.rel <- cbind(rsphere_stage4_nm.df, rsphere_4_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rsphere_stage4_nm.rel[sapply(rsphere_stage4_nm.rel, is.character)] <- lapply(rsphere_stage4_nm.rel[sapply(rsphere_stage4_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Early Reproductive Indicator Taxa Data Frame
# Subset indicator taxa data frame to only early reproductive indicators
rsphere_early.reprod_indic <- indic_rsphere_stage_aggre_100.tdy[indic_rsphere_stage_aggre_100.tdy$stage_aggregate == "early_reproductive", ]
dim(rsphere_early.reprod_indic) # 100 early reproductive indicator taxa
# Check if all early reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rsphere_stage4_nm.rel <- rownames_to_column(rsphere_stage4_nm.rel, var = "asvID")
sum(rsphere_early.reprod_indic$asvID %in% rsphere_stage4_nm.rel$asvID) == nrow(rsphere_early.reprod_indic) # TRUE → all early reproductive indicator ASVs from the "rsphere_early.reprod_indic" are present in the "rsphere_stage4_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rsphere_stage4_nm.rel.indic <- merge.data.frame(rsphere_stage4_nm.rel, rsphere_early.reprod_indic, by="asvID", all.x = T)
# Check
sum(rsphere_early.reprod_indic$asvID %in% rsphere_stage4_nm.rel$asvID) # n = 100, This gives the number of ASVs from "rsphere_early.reprod_indic" that exist in "rsphere_stage4_nm.rel".
sum(!rsphere_early.reprod_indic$asvID %in% rsphere_stage4_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Early Reproductive
rsphere_stage4_nm.rel.indic2 <- rsphere_stage4_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rsphere.early.reprod.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rsphere_stage4_nm.rel.indic$asvID, rsphere_stage4_nm.rel.indic2$asvID)
# Check if all early reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rsphere_early.reprod_indic$asvID %in% rsphere_stage4_nm.rel.indic2$asvID)
sum(rsphere_early.reprod_indic$asvID %in% 
      rsphere_stage4_nm.rel.indic2$asvID) == nrow(rsphere_early.reprod_indic) #TRUE
# Adding stage indicator column (TRUE or FALSE)
rsphere_stage4_nm.rel.indic2 <- rsphere_stage4_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Early Reproductive taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rsphere_stage4_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Early_Reproductive_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rsphere.above.stage4 <- subset(rsphere_stage4_nm.rel.indic2, above_CI == TRUE)
dim(rsphere.above.stage4) # 342 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Early reproductive
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = early_reproductive ; stage_indicator = TRUE
rsphere_stage4_nm.rel.indic2 <- column_to_rownames(rsphere_stage4_nm.rel.indic2, var = "asvID")
rsphere.indicsp.sncm.overlap.early.reprod <- subset(rsphere_stage4_nm.rel.indic2,
                                                    above_CI == TRUE &
                                                      stage_aggregate == "early_reproductive" &
                                                      stage_indicator == TRUE)
dim(rsphere.indicsp.sncm.overlap.early.reprod) # 14 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rsphere.indicsp.sncm.overlap.early.reprod, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizosphere_selected_indicator_early.reproductive.csv")

### 5.) Late Reproductive
# a) Calculate log Mean Relative Abundance and Occupancy 
# Take meta data from the late reproductive phyloseq object
rsphere_5_metadata <- data.frame(sample_data(rsphere_5_ps))
# Compute mean relative abundance
rsphere_5_asv <- otu_table(rsphere_5_ps)
rsphere_5_asv_rel <- apply(decostand(rsphere_5_asv, method="total", MARGIN=2),1, mean)
rsphere_5_asv_rel <- data.frame(rsphere_5_asv_rel)
# Compute log10(mean relative abundance)
rsphere_5_asv_rel$log_mean <- log10(rsphere_5_asv_rel$rsphere_5_asv_rel)
rsphere_5_asv_rel$asvID <- rownames(rsphere_5_asv_rel)
colnames(rsphere_5_asv_rel)[1] <- "mean_relabund"
# Compute Occupancy
rsphere_5_asv[rsphere_5_asv > 0] <- 1
rsphere_5_asv_rel$occupancy <- rowSums(rsphere_5_asv)/ncol(rsphere_5_asv)
dim(rsphere_5_asv_rel) # 19549
# b) Fit Sloan Neutral Community Model (SNCM)
# Make input data for the Sloan Neutral Model
rsphere.spp_5 <-data.frame(t(otu_table(rsphere_5_ps)), check.names = F)
rsphere.taxon_5 <- data.frame(tax_table(rsphere_5_ps))
# Fit the Sloan Neutral Community Model (SNCM) on late reproductive
set.seed(13)
rsphere_stage5_nm = fit_sncm(rsphere.spp_5, pool=NULL, rsphere.taxon_5)
rsphere_stage5_nm$fitstats
# Make the data frame of the SNCM output
rsphere_stage5_nm.df <- as.data.frame(rsphere_stage5_nm)
arrange(rsphere_stage5_nm.df, desc(predictions.freq))
rsphere_stage5_nm.df$predictions.fit_class <- as.factor(rsphere_stage5_nm.df$predictions.fit_class)
rsphere_stage5_nm.df$predictions.fit_class 
# Indicate taxa that are above CI as TRUE
rsphere_stage5_nm.df$above_CI <- rsphere_stage5_nm.df$predictions.freq > rsphere_stage5_nm.df$predictions.pred.upr
dim(rsphere_stage5_nm.df) # 19549 27
# Calculate the percentage of above and below the CI of the model
rsphere.stage5.above.pred = sum(rsphere_stage5_nm.df$predictions.freq > (rsphere_stage5_nm.df$predictions.pred.upr), na.rm=TRUE)/nrow(rsphere_stage5_nm.df)  # fraction of OTUs above prediction # 1.46 %
rsphere.stage5.below.pred = sum(rsphere_stage5_nm.df$predictions.freq < (rsphere_stage5_nm.df$predictions.pred.lwr), na.rm=TRUE)/nrow(rsphere_stage5_nm.df)  # fraction of OTUs below prediction # 1.37 %
100 - (rsphere.stage5.above.pred*100) - (rsphere.stage5.below.pred*100)# 97.17%
# c) Combine the Sloan Neutral Community Model Data Frame with Relative Abundance and Occupancy Data Frame
# Check the row names order before cbind both data frames
all(rownames(rsphere_stage5_nm.df) == rownames(rsphere_5_asv_rel)) # FALSE --> need to be the same order
setequal(rownames(rsphere_stage5_nm.df), rownames(rsphere_5_asv_rel)) # not necessarily the same order
all(rownames(rsphere_stage5_nm.df) %in% rownames(rsphere_5_asv_rel))  # Should be TRUE
all(rownames(rsphere_5_asv_rel) %in% rownames(rsphere_stage5_nm.df)) # Should be TRUE
# Re-order the row names of the SNCM output data frame to match the row names of the mean relative abundance and occupancy data frame
rsphere_stage5_nm.df <- rsphere_stage5_nm.df[rownames(rsphere_5_asv_rel), , drop = FALSE] # re-order
rsphere_stage5_nm.rel <- cbind(rsphere_stage5_nm.df, rsphere_5_asv_rel[, c("mean_relabund", "log_mean", "occupancy")]) # "stage_aggregate", "stage_indicator")]) # add columns from df2 to df1
rsphere_stage5_nm.rel[sapply(rsphere_stage5_nm.rel, is.character)] <- lapply(rsphere_stage5_nm.rel[sapply(rsphere_stage5_nm.rel, is.character)], as.factor)
# d) Merge the combined SNCM, Relative Abundance, and Occupancy Data Frame with the Late Reproductive Indicator Taxa Data Frame
# Subset indicator taxa data frame to only late reproductive indicators
rsphere_late.reprod_indic <- indic_rsphere_stage_aggre_100.tdy[indic_rsphere_stage_aggre_100.tdy$stage_aggregate == "late_reproductive", ]
dim(rsphere_late.reprod_indic) # 307 late reproductive indicator taxa
# Check if all late reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
rsphere_stage5_nm.rel <- rownames_to_column(rsphere_stage5_nm.rel, var = "asvID")
sum(rsphere_late.reprod_indic$asvID %in% rsphere_stage5_nm.rel$asvID) == nrow(rsphere_late.reprod_indic) # TRUE → all late reproductive indicator ASVs from the "rsphere_late.reprod_indic" are present in the "rsphere_stage5_nm.rel"
# Merge the mean relative abundance and occupancy data frame with the indicator taxa data frame from the IndVal analysis 
rsphere_stage5_nm.rel.indic <- merge.data.frame(rsphere_stage5_nm.rel, rsphere_late.reprod_indic, by="asvID", all.x = T)
# Check
sum(rsphere_late.reprod_indic$asvID %in% rsphere_stage5_nm.rel$asvID) # n = 307, This gives the number of ASVs from "rplane_late.reprod_indic" that exist in "rsphere_stage5_nm.rel".
sum(!rsphere_late.reprod_indic$asvID %in% rsphere_stage5_nm.rel$asvID) # 0
# add IndVal.Stat for non indicator ASVs in Late Reproductive
rsphere_stage5_nm.rel.indic2 <- rsphere_stage5_nm.rel.indic %>%
  select(-IndVal.Stat, -IndVal.p_value) %>%
  left_join(rsphere.late.reprod.indval %>%
              select(asvID, IndVal.Stat, IndVal.p_value), by = "asvID") %>%
  relocate(IndVal.Stat, IndVal.p_value, .before = IndVal.signif)
setequal(rsphere_stage5_nm.rel.indic$asvID, rsphere_stage5_nm.rel.indic2$asvID)
# Check if all late reproductive indicator ASVs are present in the SNCM relative abundance and occupancy data frame
all(rsphere_late.reprod_indic$asvID %in% rsphere_stage5_nm.rel.indic2$asvID)
sum(rsphere_late.reprod_indic$asvID %in% 
      rsphere_stage5_nm.rel.indic2$asvID) == nrow(rsphere_late.reprod_indic) #TRUE
# Adding stage indicator column (TRUE or FALSE)
rsphere_stage5_nm.rel.indic2 <- rsphere_stage5_nm.rel.indic2 %>% 
  mutate(stage_indicator = if_else(is.na(stage_aggregate), FALSE, TRUE))
# Save the complete Late Reproductive taxa that contain these information: Indicators and their indicator values, Mean Relative Abundance, Occupancy, and Neutral Model prediction
#write.csv(rsphere_stage5_nm.rel.indic2, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Complete_Data_by_Stage/Rhizosphere_Late_Reproductive_complete_SNCM_Indic.csv")
# Subset data with ASVs that are above CI = deterministic = selected by plant
rsphere.above.stage5 <- subset(rsphere_stage5_nm.rel.indic2, above_CI == TRUE)
dim(rsphere.above.stage5) # 286 taxa are above CI
# e) Identify overlapping/shared taxa between indicator analysis and SNCM (above CI) in Late reproductive
# Subset the data that meet these conditions : above_CI = TRUE ; stage_aggregate = late_reproductive ; stage_indicator = TRUE
rsphere_stage5_nm.rel.indic2 <- column_to_rownames(rsphere_stage5_nm.rel.indic2, var = "asvID")
rsphere.indicsp.sncm.overlap.late.reprod <- subset(rsphere_stage5_nm.rel.indic2,
                                                   above_CI == TRUE &
                                                     stage_aggregate == "late_reproductive" &
                                                     stage_indicator == TRUE)
dim(rsphere.indicsp.sncm.overlap.late.reprod) # 38 Host-Selected Taxa
# Save the Host-Selected Taxa
#write.csv(rsphere.indicsp.sncm.overlap.late.reprod, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Sloan_Indicator_selected/rhizosphere_selected_indicator_late.reproductive.csv")


################################################################################################
### Combine All Host-Selected Taxa from All Growth Stages and Add Updated SILVA 144 Taxonomy ###
################################################################################################

# Read the updated Uniform Taxonomy SILVA 144
new.tax.unif <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/uniform_taxonomy_edit.csv", row.names = 1)
new.tax.unif <- rownames_to_column(new.tax.unif, var = "asvID")
# Read the updated Weighted Taxonomy SILVA 144 (Weighted = based on specific habitat/sampling location, here I used the weighted taxonomy of Plant-Rhizosphere Soil = https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/V3V4-341f-806r/weighted/plant-rhizosphere)
new.tax.weight <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/weight_taxonomy_edit.csv", row.names = 1)
new.tax.weight <- rownames_to_column(new.tax.weight, var = "asvID")


# 1. RHIZOPLANE


rplane_INdSN.sharedASVs_stage1 <- rp.indicsp.sncm.overlap.seedling
rplane_INdSN.sharedASVs_stage2 <- rp.indicsp.sncm.overlap.early.veget
rplane_INdSN.sharedASVs_stage3 <- rp.indicsp.sncm.overlap.late.veget
rplane_INdSN.sharedASVs_stage4 <- rp.indicsp.sncm.overlap.early.reprod
rplane_INdSN.sharedASVs_stage5 <- rp.indicsp.sncm.overlap.late.reprod
# Adding column stage with proper stage name
rplane_INdSN.sharedASVs_stage1$Stage <- "Seedling"
rplane_INdSN.sharedASVs_stage2$Stage <- "Early Vegetative"
rplane_INdSN.sharedASVs_stage3$Stage <- "Late Vegetative"
rplane_INdSN.sharedASVs_stage4$Stage <- "Early Reproductive"
rplane_INdSN.sharedASVs_stage5$Stage <- "Late Reproductive"
# Combine all host-selected taxa data frames
rp.above_Ind <- bind_rows(rplane_INdSN.sharedASVs_stage1,rplane_INdSN.sharedASVs_stage2,rplane_INdSN.sharedASVs_stage3,rplane_INdSN.sharedASVs_stage4,rplane_INdSN.sharedASVs_stage5)
dim(rp.above_Ind) # 159 ASVs
View(rp.above_Ind)
# Re-name the old taxonomy SILVA 138 column name
rp.above_Ind.edtax <- rp.above_Ind %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rp.above_Ind.edtax <- rownames_to_column(rp.above_Ind.edtax, var = "asvID")
rp.above_Ind.edtax <- rp.above_Ind.edtax %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new SILVA 144 taxonomy columns
rp.above_Ind.edtax <- rp.above_Ind.edtax %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Edit the taxonomy columns
rp.above_Ind.edtax <- rp.above_Ind.edtax %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__")) %>%
  mutate(Weight.Class = str_remove(Weight.Class, "^c__")) %>%
  mutate(Weight.Order = str_remove(Weight.Order, "^o__")) %>%
  mutate(Weight.Family = str_remove(Weight.Family, "^f__")) %>%
  mutate(Weight.Genus = str_remove(Weight.Genus, "^g__")) %>%
  mutate(Weight.Species = str_remove(Weight.Species, "^s__")) 
# Make a new column with combination between Weight.Family and ASV ID for plotting
rp.above_Ind.edtax <- rp.above_Ind.edtax %>%
  mutate(Weight.Family_asvID = paste0(case_when(
    !is.na(Weight.Family) & Weight.Family != "" & Weight.Family != "--" & Weight.Family != "Incertae_Sedis" ~ Weight.Family,
    !is.na(Weight.Order) & Weight.Order != "" & Weight.Order != "--" & Weight.Order != "Incertae_Sedis" ~ Weight.Order,
    !is.na(Weight.Class) & Weight.Class != "" & Weight.Class != "--" & Weight.Class != "Incertae_Sedis" ~ Weight.Class,
    !is.na(Weight.Phylum) & Weight.Phylum != "" & Weight.Phylum != "--" & Weight.Phylum != "Incertae_Sedis" ~ Weight.Phylum, 
    TRUE ~ NA_character_), "-", substr(asvID, 1, 4))) 
# Save in the computer
#write.csv(rp.above_Ind.edtax, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Rhizoplane_Neutral_Indicator_SharedASVs.csv")


# 2. RHIZOSPHERE


rsphere_INdSN.sharedASVs_stage1 <- rsphere.indicsp.sncm.overlap.seedling
rsphere_INdSN.sharedASVs_stage2 <- rsphere.indicsp.sncm.overlap.early.veget
rsphere_INdSN.sharedASVs_stage3 <- rsphere.indicsp.sncm.overlap.late.veget
rsphere_INdSN.sharedASVs_stage4 <- rsphere.indicsp.sncm.overlap.early.reprod
rsphere_INdSN.sharedASVs_stage5 <- rsphere.indicsp.sncm.overlap.late.reprod
# Adding column stage with proper stage name
rsphere_INdSN.sharedASVs_stage1$Stage <- "Seedling"
rsphere_INdSN.sharedASVs_stage2$Stage <- "Early Vegetative"
rsphere_INdSN.sharedASVs_stage3$Stage <- "Late Vegetative"
rsphere_INdSN.sharedASVs_stage4$Stage <- "Early Reproductive"
rsphere_INdSN.sharedASVs_stage5$Stage <- "Late Reproductive"
# Combine all host-selected taxa data frames
rs.above_Ind <- bind_rows(rsphere_INdSN.sharedASVs_stage1, rsphere_INdSN.sharedASVs_stage2, 
                          rsphere_INdSN.sharedASVs_stage3, rsphere_INdSN.sharedASVs_stage4,
                          rsphere_INdSN.sharedASVs_stage5)
dim(rs.above_Ind) # 131 ASVs
# Re-name the old taxonomy SILVA 138 column name
rs.above_Ind.edtax <- rs.above_Ind %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rs.above_Ind.edtax <- rownames_to_column(rs.above_Ind.edtax, var = "asvID")
rs.above_Ind.edtax <- rs.above_Ind.edtax %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new SILVA 144 taxonomy columns
rs.above_Ind.edtax <- rs.above_Ind.edtax %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Edit the taxonomy columns
rs.above_Ind.edtax <- rs.above_Ind.edtax %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__")) %>%
  mutate(Weight.Class = str_remove(Weight.Class, "^c__")) %>%
  mutate(Weight.Order = str_remove(Weight.Order, "^o__")) %>%
  mutate(Weight.Family = str_remove(Weight.Family, "^f__")) %>%
  mutate(Weight.Genus = str_remove(Weight.Genus, "^g__")) %>%
  mutate(Weight.Species = str_remove(Weight.Species, "^s__")) 
# Make a new column with combination between Weight.Family and ASV ID for plotting
rs.above_Ind.edtax <- rs.above_Ind.edtax %>%
  mutate(Weight.Family_asvID = paste0(case_when(
    !is.na(Weight.Family) & Weight.Family != "" & Weight.Family != "--" & Weight.Family != "Incertae_Sedis" ~ Weight.Family,
    !is.na(Weight.Order) & Weight.Order != "" & Weight.Order != "--" & Weight.Order != "Incertae_Sedis" ~ Weight.Order,
    !is.na(Weight.Class) & Weight.Class != "" & Weight.Class != "--" & Weight.Class != "Incertae_Sedis" ~ Weight.Class,
    !is.na(Weight.Phylum) & Weight.Phylum != "" & Weight.Phylum != "--" & Weight.Phylum != "Incertae_Sedis" ~ Weight.Phylum, 
    TRUE ~ NA_character_), "-", substr(asvID, 1, 4))) 
# Save in the computer
#write.csv(rs.above_Ind.edtax, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/Rhizosphere_Neutral_Indicator_SharedASVs.csv")



#########################################
### Tidy Up All Data Set for Plotting ###
#########################################

### 1. Tidy up the Growth Stage Complete Dataset 

### I. RHIZOPLANE

# 1.) Load the complete Seedling data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rplane_stage1_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizoplane_Seedling_complete_SNCM_Indic.csv", row.names = 1)
head(rplane_stage1_nm.rel.indic)
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
rplane_stage1_nm.rel.indic.df$Stage <- "Seedling"
rplane_stage1_nm.rel.indic.df$stage_aggregate <- as.factor(rplane_stage1_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_stage1_nm.rel.indic.df <- rplane_stage1_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 2.) Load Early Vegetative Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rplane_stage2_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizoplane_Early_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
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
rplane_stage2_nm.rel.indic.df$Stage <- "Early Vegetative"
rplane_stage2_nm.rel.indic.df$stage_aggregate <- as.factor(rplane_stage2_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_stage2_nm.rel.indic.df <- rplane_stage2_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 3.) Load Late Vegetative Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rplane_stage3_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizoplane_Late_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
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
rplane_stage3_nm.rel.indic.df$Stage <- "Late Vegetative"
rplane_stage3_nm.rel.indic.df$stage_aggregate <- as.factor(rplane_stage3_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_stage3_nm.rel.indic.df <- rplane_stage3_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 4.) Load Early Reproductive Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rplane_stage4_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizoplane_Early_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
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
rplane_stage4_nm.rel.indic.df$Stage <- "Early Reproductive"
rplane_stage4_nm.rel.indic.df$stage_aggregate <- as.factor(rplane_stage4_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_stage4_nm.rel.indic.df <- rplane_stage4_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 5.) Load Late Reproductive Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rplane_stage5_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizoplane_Late_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
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
rplane_stage5_nm.rel.indic.df$Stage <- "Late Reproductive"
rplane_stage5_nm.rel.indic.df$stage_aggregate <- as.factor(rplane_stage5_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_stage5_nm.rel.indic.df <- rplane_stage5_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)


### II. RHIZOSPHERE

# 1.) Load Seedling Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rsphere_stage1_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizosphere_Seedling_complete_SNCM_Indic.csv", row.names = 1)
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
rsphere_stage1_nm.rel.indic.df$Stage <- "Seedling"
rsphere_stage1_nm.rel.indic.df$stage_aggregate <- as.factor(rsphere_stage1_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_stage1_nm.rel.indic.df <- rsphere_stage1_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 2.) Load Early Vegetative Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rsphere_stage2_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizosphere_Early_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
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
rsphere_stage2_nm.rel.indic.df$Stage <- "Early Vegetative"
rsphere_stage2_nm.rel.indic.df$stage_aggregate <- as.factor(rsphere_stage2_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_stage2_nm.rel.indic.df <- rsphere_stage2_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 3.) Load Late Vegetative Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rsphere_stage3_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizosphere_Late_Vegetative_complete_SNCM_Indic.csv", row.names = 1)
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
rsphere_stage3_nm.rel.indic.df$Stage <- "Late Vegetative"
rsphere_stage3_nm.rel.indic.df$stage_aggregate <- as.factor(rsphere_stage3_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_stage3_nm.rel.indic.df <- rsphere_stage3_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 4.) Load Early Reproductive Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rsphere_stage4_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage/Rhizosphere_Early_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
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
rsphere_stage4_nm.rel.indic.df$Stage <- "Early Reproductive"
rsphere_stage4_nm.rel.indic.df$stage_aggregate <- as.factor(rsphere_stage4_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_stage4_nm.rel.indic.df <- rsphere_stage4_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)

# 5.) Load Late Reproductive Data set that already contain Indicator values and Neutral Model prediction (as well as the mean relative abundance and occupancy)
rsphere_stage5_nm.rel.indic <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Complete_Data_by_Stage//Rhizosphere_Late_Reproductive_complete_SNCM_Indic.csv", row.names = 1)
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
rsphere_stage5_nm.rel.indic.df$Stage <- "Late Reproductive"
rsphere_stage5_nm.rel.indic.df$stage_aggregate <- as.factor(rsphere_stage5_nm.rel.indic.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_stage5_nm.rel.indic.df <- rsphere_stage5_nm.rel.indic.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)



### 2. Build a complete dataset that includes indicator taxa for each growth stage, as well as the same taxa that are present (or absent) in other stages where they are not indicators.


### I. RHIZOPLANE

# Make asv ID as row names in all growth stage dataset
rplane_stage1_nm.rel.indic.df <- column_to_rownames(rplane_stage1_nm.rel.indic.df, var = "asvID")
rplane_stage2_nm.rel.indic.df <- column_to_rownames(rplane_stage2_nm.rel.indic.df, var = "asvID")
rplane_stage3_nm.rel.indic.df <- column_to_rownames(rplane_stage3_nm.rel.indic.df, var = "asvID")
rplane_stage4_nm.rel.indic.df <- column_to_rownames(rplane_stage4_nm.rel.indic.df, var = "asvID")
rplane_stage5_nm.rel.indic.df <- column_to_rownames(rplane_stage5_nm.rel.indic.df, var = "asvID")

# 1.) Seedling
# Load the indicator taxa which fall in the above prediction in the Neutral model
rplane_INdSN.sharedASVs_stage1 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizoplane_selected_indicator_seedling.csv")
colnames(rplane_INdSN.sharedASVs_stage1) # 30
# Re-name the old taxonomy SILVA 138 column name
rplane_INdSN.sharedASVs_stage1 <- rplane_INdSN.sharedASVs_stage1 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_INdSN.sharedASVs_stage1 <- rplane_INdSN.sharedASVs_stage1 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_INdSN.sharedASVs_stage1 <- rplane_INdSN.sharedASVs_stage1 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_INdSN.sharedASVs_stage1.df <- rplane_INdSN.sharedASVs_stage1 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_INdSN.sharedASVs_stage1.df$point_class <- ifelse(rplane_INdSN.sharedASVs_stage1.df$above_CI == TRUE,
                                                        "Above",ifelse(rplane_INdSN.sharedASVs_stage1.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_INdSN.sharedASVs_stage1.df$point_class <- factor(rplane_INdSN.sharedASVs_stage1.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rplane_INdSN.sharedASVs_stage1.df$Stage <- "Seedling"
rplane_INdSN.sharedASVs_stage1.df$stage_aggregate <- as.factor(rplane_INdSN.sharedASVs_stage1.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_INdSN.sharedASVs_stage1.df <- rplane_INdSN.sharedASVs_stage1.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rplane_INdSN.sharedASVs_stage1.df <- column_to_rownames(rplane_INdSN.sharedASVs_stage1.df, var = "asvID")
# Subset those Seedling selected ASVs from the other stages:
# - Early vegetative
rp.early.veget_for.seedling.all <- rplane_stage2_nm.rel.indic.df[rownames(rplane_stage2_nm.rel.indic.df) %in%
                                                                   rownames(rplane_INdSN.sharedASVs_stage1.df), ]
rp.early.veget_for.seedling.all <- rownames_to_column(rp.early.veget_for.seedling.all,var = "asvID") # 28
# - Late vegetative
rp.late.veget_for.seedling.all <- rplane_stage3_nm.rel.indic.df[rownames(rplane_stage3_nm.rel.indic.df) %in%
                                                                  rownames(rplane_INdSN.sharedASVs_stage1.df), ]
rp.late.veget_for.seedling.all <- rownames_to_column(rp.late.veget_for.seedling.all,var = "asvID") # 30
# - Early reproductive
rp.early.reprod_for.seedling.all <- rplane_stage4_nm.rel.indic.df[rownames(rplane_stage4_nm.rel.indic.df) %in%
                                                                    rownames(rplane_INdSN.sharedASVs_stage1.df), ]
rp.early.reprod_for.seedling.all <- rownames_to_column(rp.early.reprod_for.seedling.all,var = "asvID") # 28
# - Late reproductive
rp.late.reprod_for.seedling.all <- rplane_stage5_nm.rel.indic.df[rownames(rplane_stage5_nm.rel.indic.df) %in%
                                                                   rownames(rplane_INdSN.sharedASVs_stage1.df), ]
rp.late.reprod_for.seedling.all <- rownames_to_column(rp.late.reprod_for.seedling.all,var = "asvID") # 28
# Combine those selected seedling ASVs with the same ASVs from other stages
rplane_INdSN.sharedASVs_stage1.df <- rownames_to_column(rplane_INdSN.sharedASVs_stage1.df, var = "asvID")
rp.seedling.selected.all <- bind_rows(rplane_INdSN.sharedASVs_stage1.df, rp.early.veget_for.seedling.all, rp.late.veget_for.seedling.all,
                                      rp.early.reprod_for.seedling.all,rp.late.reprod_for.seedling.all)

# 2.) Early Vegetative
# Load the indicator taxa which fall in the above prediction in the Neutral model
rplane_INdSN.sharedASVs_stage2 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizoplane_selected_indicator_early.vegetative.csv")
dim(rplane_INdSN.sharedASVs_stage2) # 25
# Re-name the old taxonomy SILVA 138 column name
rplane_INdSN.sharedASVs_stage2 <- rplane_INdSN.sharedASVs_stage2 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_INdSN.sharedASVs_stage2 <- rplane_INdSN.sharedASVs_stage2 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_INdSN.sharedASVs_stage2 <- rplane_INdSN.sharedASVs_stage2 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_INdSN.sharedASVs_stage2.df <- rplane_INdSN.sharedASVs_stage2 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_INdSN.sharedASVs_stage2.df$point_class <- ifelse(rplane_INdSN.sharedASVs_stage2.df$above_CI == TRUE,
                                                        "Above",ifelse(rplane_INdSN.sharedASVs_stage2.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_INdSN.sharedASVs_stage2.df$point_class <- factor(rplane_INdSN.sharedASVs_stage2.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rplane_INdSN.sharedASVs_stage2.df$Stage <- "Early Vegetative"
rplane_INdSN.sharedASVs_stage2.df$stage_aggregate <- as.factor(rplane_INdSN.sharedASVs_stage2.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_INdSN.sharedASVs_stage2.df <- rplane_INdSN.sharedASVs_stage2.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)
rplane_INdSN.sharedASVs_stage2.df <- column_to_rownames(rplane_INdSN.sharedASVs_stage2.df, var = "asvID")
# Subset those Early Vegetative selected ASVs from the other stages:
# - Seedling
rp.seedling_for.early.veget.all <- rplane_stage1_nm.rel.indic.df[rownames(rplane_stage1_nm.rel.indic.df) %in%
                                                                   rownames(rplane_INdSN.sharedASVs_stage2.df), ]
rp.seedling_for.early.veget.all <- rownames_to_column(rp.seedling_for.early.veget.all,var = "asvID") # 22
# - Late vegetative
rp.late.veget_for.early.veget.all <- rplane_stage3_nm.rel.indic.df[rownames(rplane_stage3_nm.rel.indic.df) %in%
                                                                     rownames(rplane_INdSN.sharedASVs_stage2.df), ]
rp.late.veget_for.early.veget.all <- rownames_to_column(rp.late.veget_for.early.veget.all,var = "asvID") # 24
# - Early reproductive
rp.early.reprod_for.early.veget.all <- rplane_stage4_nm.rel.indic.df[rownames(rplane_stage4_nm.rel.indic.df) %in%
                                                                       rownames(rplane_INdSN.sharedASVs_stage2.df), ]
rp.early.reprod_for.early.veget.all <- rownames_to_column(rp.early.reprod_for.early.veget.all,var = "asvID") # 24
# - Late reproductive
rp.late.reprod_for.early.veget.all <- rplane_stage5_nm.rel.indic.df[rownames(rplane_stage5_nm.rel.indic.df) %in%
                                                                      rownames(rplane_INdSN.sharedASVs_stage2.df), ]
rp.late.reprod_for.early.veget.all <- rownames_to_column(rp.late.reprod_for.early.veget.all,var = "asvID") # 19
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rplane_INdSN.sharedASVs_stage2.df <- rownames_to_column(rplane_INdSN.sharedASVs_stage2.df, var = "asvID")
rp.early.veget.selected.all <- bind_rows(rp.seedling_for.early.veget.all, rplane_INdSN.sharedASVs_stage2.df, rp.late.veget_for.early.veget.all,
                                         rp.early.reprod_for.early.veget.all, rp.late.reprod_for.early.veget.all)

# 3.) Late Vegetative
# Load the indicator taxa which fall in the above prediction in the Neutral model
rplane_INdSN.sharedASVs_stage3 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizoplane_selected_indicator_late.vegetative.csv")
dim(rplane_INdSN.sharedASVs_stage3) # 54
# Re-name the old taxonomy SILVA 138 column name
rplane_INdSN.sharedASVs_stage3 <- rplane_INdSN.sharedASVs_stage3 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_INdSN.sharedASVs_stage3 <- rplane_INdSN.sharedASVs_stage3 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_INdSN.sharedASVs_stage3 <- rplane_INdSN.sharedASVs_stage3 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_INdSN.sharedASVs_stage3.df <- rplane_INdSN.sharedASVs_stage3 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_INdSN.sharedASVs_stage3.df$point_class <- ifelse(rplane_INdSN.sharedASVs_stage3.df$above_CI == TRUE,
                                                        "Above",ifelse(rplane_INdSN.sharedASVs_stage3.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_INdSN.sharedASVs_stage3.df$point_class <- factor(rplane_INdSN.sharedASVs_stage3.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rplane_INdSN.sharedASVs_stage3.df$Stage <- "Late Vegetative"
rplane_INdSN.sharedASVs_stage3.df$stage_aggregate <- as.factor(rplane_INdSN.sharedASVs_stage3.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_INdSN.sharedASVs_stage3.df <- rplane_INdSN.sharedASVs_stage3.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)
rplane_INdSN.sharedASVs_stage3.df <- column_to_rownames(rplane_INdSN.sharedASVs_stage3.df, var = "asvID")
# Subset those Late Vegetative selected ASVs from the other stages:
# - Seedling
rp.seedling_for.late.veget.all <- rplane_stage1_nm.rel.indic.df[rownames(rplane_stage1_nm.rel.indic.df) %in%
                                                                  rownames(rplane_INdSN.sharedASVs_stage3.df), ]
rp.seedling_for.late.veget.all <- rownames_to_column(rp.seedling_for.late.veget.all,var = "asvID") # 51
# - Early vegetative
rp.early.veget_for.late.veget.all <- rplane_stage2_nm.rel.indic.df[rownames(rplane_stage2_nm.rel.indic.df) %in%
                                                                     rownames(rplane_INdSN.sharedASVs_stage3.df), ]
rp.early.veget_for.late.veget.all <- rownames_to_column(rp.early.veget_for.late.veget.all,var = "asvID") # 49
# - Early reproductive
rp.early.reprod_for.late.veget.all <- rplane_stage4_nm.rel.indic.df[rownames(rplane_stage4_nm.rel.indic.df) %in%
                                                                      rownames(rplane_INdSN.sharedASVs_stage3.df), ]
rp.early.reprod_for.late.veget.all <- rownames_to_column(rp.early.reprod_for.late.veget.all,var = "asvID") # 53
# - Late reproductive
rp.late.reprod_for.late.veget.all <- rplane_stage5_nm.rel.indic.df[rownames(rplane_stage5_nm.rel.indic.df) %in%
                                                                     rownames(rplane_INdSN.sharedASVs_stage3.df), ]
rp.late.reprod_for.late.veget.all <- rownames_to_column(rp.late.reprod_for.late.veget.all,var = "asvID") # 50
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rplane_INdSN.sharedASVs_stage3.df <- rownames_to_column(rplane_INdSN.sharedASVs_stage3.df, var = "asvID")
rp.late.veget.selected.all <- bind_rows(rp.seedling_for.late.veget.all, rp.early.veget_for.late.veget.all, rplane_INdSN.sharedASVs_stage3.df,
                                        rp.early.reprod_for.late.veget.all, rp.late.reprod_for.late.veget.all)

# 4.) Early Reproductive
# Load the indicator taxa which fall in the above prediction in the Neutral model
rplane_INdSN.sharedASVs_stage4 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizoplane_selected_indicator_early.reproductive.csv")
dim(rplane_INdSN.sharedASVs_stage4) # 26
# Re-name the old taxonomy SILVA 138 column name
rplane_INdSN.sharedASVs_stage4 <- rplane_INdSN.sharedASVs_stage4 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_INdSN.sharedASVs_stage4 <- rplane_INdSN.sharedASVs_stage4 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_INdSN.sharedASVs_stage4 <- rplane_INdSN.sharedASVs_stage4 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_INdSN.sharedASVs_stage4.df <- rplane_INdSN.sharedASVs_stage4 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_INdSN.sharedASVs_stage4.df$point_class <- ifelse(rplane_INdSN.sharedASVs_stage4.df$above_CI == TRUE,
                                                        "Above",ifelse(rplane_INdSN.sharedASVs_stage4.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_INdSN.sharedASVs_stage4.df$point_class <- factor(rplane_INdSN.sharedASVs_stage4.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rplane_INdSN.sharedASVs_stage4.df$Stage <- "Early Reproductive"
rplane_INdSN.sharedASVs_stage4.df$stage_aggregate <- as.factor(rplane_INdSN.sharedASVs_stage4.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_INdSN.sharedASVs_stage4.df <- rplane_INdSN.sharedASVs_stage4.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)
rplane_INdSN.sharedASVs_stage4.df <- column_to_rownames(rplane_INdSN.sharedASVs_stage4.df, var = "asvID")
# Subset those Early Reproductive selected ASVs from the other stages:
# - Seedling
rp.seedling_for.early.reprod.all <- rplane_stage1_nm.rel.indic.df[rownames(rplane_stage1_nm.rel.indic.df) %in%
                                                                    rownames(rplane_INdSN.sharedASVs_stage4.df), ]
rp.seedling_for.early.reprod.all <- rownames_to_column(rp.seedling_for.early.reprod.all,var = "asvID") # 26
# - Early vegetative
rp.early.veget_for.early.reprod.all <- rplane_stage2_nm.rel.indic.df[rownames(rplane_stage2_nm.rel.indic.df) %in%
                                                                       rownames(rplane_INdSN.sharedASVs_stage4.df), ]
rp.early.veget_for.early.reprod.all <- rownames_to_column(rp.early.veget_for.early.reprod.all,var = "asvID") # 24
# - Late vegetative
rp.late.veget_for.early.reprod.all <- rplane_stage3_nm.rel.indic.df[rownames(rplane_stage3_nm.rel.indic.df) %in%
                                                                      rownames(rplane_INdSN.sharedASVs_stage4.df), ]
rp.late.veget_for.early.reprod.all <- rownames_to_column(rp.late.veget_for.early.reprod.all,var = "asvID") # 26
# - Late reproductive
rp.late.reprod_for.early.reprod.all <- rplane_stage5_nm.rel.indic.df[rownames(rplane_stage5_nm.rel.indic.df) %in%
                                                                       rownames(rplane_INdSN.sharedASVs_stage4.df), ]
rp.late.reprod_for.early.reprod.all <- rownames_to_column(rp.late.reprod_for.early.reprod.all,var = "asvID") # 23
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rplane_INdSN.sharedASVs_stage4.df <- rownames_to_column(rplane_INdSN.sharedASVs_stage4.df, var = "asvID")
rp.early.reprod.selected.all <- bind_rows(rp.seedling_for.early.reprod.all, rp.early.veget_for.early.reprod.all, rp.late.veget_for.early.reprod.all,
                                          rplane_INdSN.sharedASVs_stage4.df, rp.late.reprod_for.early.reprod.all)

# 5.) Late Reproductive
# Load the indicator taxa which fall in the above prediction in the Neutral model
rplane_INdSN.sharedASVs_stage5 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizoplane_selected_indicator_late.reproductive.csv")
dim(rplane_INdSN.sharedASVs_stage5) # 24
# Re-name the old taxonomy SILVA 138 column name
rplane_INdSN.sharedASVs_stage5 <- rplane_INdSN.sharedASVs_stage5 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rplane_INdSN.sharedASVs_stage5 <- rplane_INdSN.sharedASVs_stage5 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rplane_INdSN.sharedASVs_stage5 <- rplane_INdSN.sharedASVs_stage5 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rplane_INdSN.sharedASVs_stage5.df <- rplane_INdSN.sharedASVs_stage5 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rplane_INdSN.sharedASVs_stage5.df$point_class <- ifelse(rplane_INdSN.sharedASVs_stage5.df$above_CI == TRUE,
                                                        "Above",ifelse(rplane_INdSN.sharedASVs_stage5.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rplane_INdSN.sharedASVs_stage5.df$point_class <- factor(rplane_INdSN.sharedASVs_stage5.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rplane_INdSN.sharedASVs_stage5.df$Stage <- "Late Reproductive"
rplane_INdSN.sharedASVs_stage5.df$stage_aggregate <- as.factor(rplane_INdSN.sharedASVs_stage5.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rplane_INdSN.sharedASVs_stage5.df <- rplane_INdSN.sharedASVs_stage5.df %>%
  mutate(mean_relabund_percent = mean_relabund * 100) %>%
  relocate(mean_relabund_percent, .after = mean_relabund)
rplane_INdSN.sharedASVs_stage5.df <- column_to_rownames(rplane_INdSN.sharedASVs_stage5.df, var = "asvID")
# Subset those Late Reproductive selected ASVs from the other stages:
# - Seedling
rp.seedling_for.late.reprod.all <- rplane_stage1_nm.rel.indic.df[rownames(rplane_stage1_nm.rel.indic.df) %in%
                                                                   rownames(rplane_INdSN.sharedASVs_stage5.df), ]
rp.seedling_for.late.reprod.all <- rownames_to_column(rp.seedling_for.late.reprod.all,var = "asvID") # 22
# - Early vegetative
rp.early.veget_for.late.reprod.all <- rplane_stage2_nm.rel.indic.df[rownames(rplane_stage2_nm.rel.indic.df) %in%
                                                                      rownames(rplane_INdSN.sharedASVs_stage5.df), ]
rp.early.veget_for.late.reprod.all <- rownames_to_column(rp.early.veget_for.late.reprod.all,var = "asvID") # 19
# - Late vegetative
rp.late.veget_for.late.reprod.all <- rplane_stage3_nm.rel.indic.df[rownames(rplane_stage3_nm.rel.indic.df) %in%
                                                                     rownames(rplane_INdSN.sharedASVs_stage5.df), ]
rp.late.veget_for.late.reprod.all <- rownames_to_column(rp.late.veget_for.late.reprod.all,var = "asvID") # 18
# - Early reproductive
rp.early.reprod_for.late.reprod.all <- rplane_stage4_nm.rel.indic.df[rownames(rplane_stage4_nm.rel.indic.df) %in%
                                                                       rownames(rplane_INdSN.sharedASVs_stage5.df), ]
rp.early.reprod_for.late.reprod.all <- rownames_to_column(rp.early.reprod_for.late.reprod.all,var = "asvID") # 21
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rplane_INdSN.sharedASVs_stage5.df <- rownames_to_column(rplane_INdSN.sharedASVs_stage5.df, var = "asvID")
rp.late.reprod.selected.all <- bind_rows(rp.seedling_for.late.reprod.all, rp.early.veget_for.late.reprod.all, rp.late.veget_for.late.reprod.all,
                                         rp.early.reprod_for.late.reprod.all, rplane_INdSN.sharedASVs_stage5.df)

### Combine all complete dataset that includes indicator taxa for each growth stage, as well as the same taxa that are present (or absent) in other stages where they are not indicators.
rp.indic.above.all <- bind_rows(rp.seedling.selected.all, 
                                rp.early.veget.selected.all, 
                                rp.late.veget.selected.all,
                                rp.early.reprod.selected.all, 
                                rp.late.reprod.selected.all)
# Save the data and make edits on the local computer
#write.csv(rp.indic.above.all, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizoplane_indic_above_all.csv")
View(rp.indic.above.all)



### II. RHIZOSHERE

# Make asv ID as row names in all growth stage dataset
rsphere_stage1_nm.rel.indic.df <- column_to_rownames(rsphere_stage1_nm.rel.indic.df, var = "asvID")
rsphere_stage2_nm.rel.indic.df <- column_to_rownames(rsphere_stage2_nm.rel.indic.df, var = "asvID")
rsphere_stage3_nm.rel.indic.df <- column_to_rownames(rsphere_stage3_nm.rel.indic.df, var = "asvID")
rsphere_stage4_nm.rel.indic.df <- column_to_rownames(rsphere_stage4_nm.rel.indic.df, var = "asvID")
rsphere_stage5_nm.rel.indic.df <- column_to_rownames(rsphere_stage5_nm.rel.indic.df, var = "asvID")

# 1.) Seedling
# Load the indicator taxa which fall in the above prediction in the Neutral model
rsphere_INdSN.sharedASVs_stage1 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizosphere_selected_indicator_seedling.csv")
dim(rsphere_INdSN.sharedASVs_stage1) # 29
# Re-name the old taxonomy SILVA 138 column name
rsphere_INdSN.sharedASVs_stage1 <- rsphere_INdSN.sharedASVs_stage1 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_INdSN.sharedASVs_stage1 <- rsphere_INdSN.sharedASVs_stage1 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_INdSN.sharedASVs_stage1 <- rsphere_INdSN.sharedASVs_stage1 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_INdSN.sharedASVs_stage1.df <- rsphere_INdSN.sharedASVs_stage1 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_INdSN.sharedASVs_stage1.df$point_class <- ifelse(rsphere_INdSN.sharedASVs_stage1.df$above_CI == TRUE,
                                                        "Above",ifelse(rsphere_INdSN.sharedASVs_stage1.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_INdSN.sharedASVs_stage1.df$point_class <- factor(rsphere_INdSN.sharedASVs_stage1.df$point_class,
                                                        levels = c("Above", "Neutral", "Below"))
rsphere_INdSN.sharedASVs_stage1.df$Stage <- "Seedling"
rsphere_INdSN.sharedASVs_stage1.df$stage_aggregate <- as.factor(rsphere_INdSN.sharedASVs_stage1.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_INdSN.sharedASVs_stage1.df <- rsphere_INdSN.sharedASVs_stage1.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rsphere_INdSN.sharedASVs_stage1.df <- column_to_rownames(rsphere_INdSN.sharedASVs_stage1.df, var = "asvID")
# Subset those Seedling selected ASVs from the other stages:
# - Early vegetative
rs.early.veget_for.seedling.all <- rsphere_stage2_nm.rel.indic.df[rownames(rsphere_stage2_nm.rel.indic.df) %in%
                                                                    rownames(rsphere_INdSN.sharedASVs_stage1.df), ]
rs.early.veget_for.seedling.all <- rownames_to_column(rs.early.veget_for.seedling.all,var = "asvID")
dim(rs.early.veget_for.seedling.all) # 28
# - Late vegetative
rs.late.veget_for.seedling.all <- rsphere_stage3_nm.rel.indic.df[rownames(rsphere_stage3_nm.rel.indic.df) %in%
                                                                   rownames(rsphere_INdSN.sharedASVs_stage1.df), ]
rs.late.veget_for.seedling.all <- rownames_to_column(rs.late.veget_for.seedling.all,var = "asvID")
dim(rs.late.veget_for.seedling.all) # 29
# - Early reproductive
rs.early.reprod_for.seedling.all <- rsphere_stage4_nm.rel.indic.df[rownames(rsphere_stage4_nm.rel.indic.df) %in%
                                                                     rownames(rsphere_INdSN.sharedASVs_stage1.df), ]
rs.early.reprod_for.seedling.all <- rownames_to_column(rs.early.reprod_for.seedling.all,var = "asvID")
dim(rs.early.reprod_for.seedling.all) # 28
# - Late reproductive
rs.late.reprod_for.seedling.all <- rsphere_stage5_nm.rel.indic.df[rownames(rsphere_stage5_nm.rel.indic.df) %in%
                                                                    rownames(rsphere_INdSN.sharedASVs_stage1.df), ]
rs.late.reprod_for.seedling.all <- rownames_to_column(rs.late.reprod_for.seedling.all,var = "asvID")
dim(rs.late.reprod_for.seedling.all) # 28
# Combine those selected seedling ASVs with the same ASVs from other stages
rsphere_INdSN.sharedASVs_stage1.df <- rownames_to_column(rsphere_INdSN.sharedASVs_stage1.df, var = "asvID")
rs.seedling.selected.all <- bind_rows(rsphere_INdSN.sharedASVs_stage1.df, rs.early.veget_for.seedling.all, rs.late.veget_for.seedling.all,
                                      rs.early.reprod_for.seedling.all, rs.late.reprod_for.seedling.all)
View(rs.seedling.selected.all)    

# 2.) Early Vegetative
# Load the indicator taxa which fall in the above prediction in the Neutral model
rsphere_INdSN.sharedASVs_stage2 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizosphere_selected_indicator_early.vegetative.csv")
dim(rsphere_INdSN.sharedASVs_stage2) # 32
# Re-name the old taxonomy SILVA 138 column name
rsphere_INdSN.sharedASVs_stage2 <- rsphere_INdSN.sharedASVs_stage2 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_INdSN.sharedASVs_stage2 <- rsphere_INdSN.sharedASVs_stage2 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_INdSN.sharedASVs_stage2 <- rsphere_INdSN.sharedASVs_stage2 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_INdSN.sharedASVs_stage2.df <- rsphere_INdSN.sharedASVs_stage2 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_INdSN.sharedASVs_stage2.df$point_class <- ifelse(rsphere_INdSN.sharedASVs_stage2.df$above_CI == TRUE,
                                                         "Above",ifelse(rsphere_INdSN.sharedASVs_stage2.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_INdSN.sharedASVs_stage2.df$point_class <- factor(rsphere_INdSN.sharedASVs_stage2.df$point_class,
                                                         levels = c("Above", "Neutral", "Below"))
rsphere_INdSN.sharedASVs_stage2.df$Stage <- "Early Vegetative"
rsphere_INdSN.sharedASVs_stage2.df$stage_aggregate <- as.factor(rsphere_INdSN.sharedASVs_stage2.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_INdSN.sharedASVs_stage2.df <- rsphere_INdSN.sharedASVs_stage2.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rsphere_INdSN.sharedASVs_stage2.df <- column_to_rownames(rsphere_INdSN.sharedASVs_stage2.df, var = "asvID")
# Subset those Early Vegetative selected ASVs from the other stages:
# - Seedling
rs.seedling_for.early.veget.all <- rsphere_stage1_nm.rel.indic.df[rownames(rsphere_stage1_nm.rel.indic.df) %in%
                                                                    rownames(rsphere_INdSN.sharedASVs_stage2.df), ]
rs.seedling_for.early.veget.all <- rownames_to_column(rs.seedling_for.early.veget.all,var = "asvID")
dim(rs.seedling_for.early.veget.all) # 32
# - Late vegetative
rs.late.veget_for.early.veget.all <- rsphere_stage3_nm.rel.indic.df[rownames(rsphere_stage3_nm.rel.indic.df) %in%
                                                                      rownames(rsphere_INdSN.sharedASVs_stage2.df), ]
rs.late.veget_for.early.veget.all <- rownames_to_column(rs.late.veget_for.early.veget.all,var = "asvID")
dim(rs.late.veget_for.early.veget.all) # 31
# - Early reproductive
rs.early.reprod_for.early.veget.all <- rsphere_stage4_nm.rel.indic.df[rownames(rsphere_stage4_nm.rel.indic.df) %in%
                                                                        rownames(rsphere_INdSN.sharedASVs_stage2.df), ]
rs.early.reprod_for.early.veget.all <- rownames_to_column(rs.early.reprod_for.early.veget.all,var = "asvID")
dim(rs.early.reprod_for.early.veget.all) # 30
# - Late reproductive
rs.late.reprod_for.early.veget.all <- rsphere_stage5_nm.rel.indic.df[rownames(rsphere_stage5_nm.rel.indic.df) %in%
                                                                       rownames(rsphere_INdSN.sharedASVs_stage2.df), ]
rs.late.reprod_for.early.veget.all <- rownames_to_column(rs.late.reprod_for.early.veget.all,var = "asvID")
dim(rs.late.reprod_for.early.veget.all) # 29
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rsphere_INdSN.sharedASVs_stage2.df <- rownames_to_column(rsphere_INdSN.sharedASVs_stage2.df, var = "asvID")
rs.early.veget.selected.all <- bind_rows(rs.seedling_for.early.veget.all, rsphere_INdSN.sharedASVs_stage2.df, rs.late.veget_for.early.veget.all,
                                         rs.early.reprod_for.early.veget.all, rs.late.reprod_for.early.veget.all)
View(rs.early.veget.selected.all)

# 3.) Late Vegetative
# Load the indicator taxa which fall in the above prediction in the Neutral model
rsphere_INdSN.sharedASVs_stage3 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizosphere_selected_indicator_late.vegetative.csv")
dim(rsphere_INdSN.sharedASVs_stage3) # 18
# Re-name the old taxonomy SILVA 138 column name
rsphere_INdSN.sharedASVs_stage3 <- rsphere_INdSN.sharedASVs_stage3 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_INdSN.sharedASVs_stage3 <- rsphere_INdSN.sharedASVs_stage3 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_INdSN.sharedASVs_stage3 <- rsphere_INdSN.sharedASVs_stage3 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_INdSN.sharedASVs_stage3.df <- rsphere_INdSN.sharedASVs_stage3 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_INdSN.sharedASVs_stage3.df$point_class <- ifelse(rsphere_INdSN.sharedASVs_stage3.df$above_CI == TRUE,
                                                         "Above",ifelse(rsphere_INdSN.sharedASVs_stage3.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_INdSN.sharedASVs_stage3.df$point_class <- factor(rsphere_INdSN.sharedASVs_stage3.df$point_class,
                                                         levels = c("Above", "Neutral", "Below"))
rsphere_INdSN.sharedASVs_stage3.df$Stage <- "Late Vegetative"
rsphere_INdSN.sharedASVs_stage3.df$stage_aggregate <- as.factor(rsphere_INdSN.sharedASVs_stage3.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_INdSN.sharedASVs_stage3.df <- rsphere_INdSN.sharedASVs_stage3.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rsphere_INdSN.sharedASVs_stage3.df <- column_to_rownames(rsphere_INdSN.sharedASVs_stage3.df, var = "asvID")
# Subset those Late Vegetative selected ASVs from the other stages:
# - Seedling
rs.seedling_for.late.veget.all <- rsphere_stage1_nm.rel.indic.df[rownames(rsphere_stage1_nm.rel.indic.df) %in%
                                                                   rownames(rsphere_INdSN.sharedASVs_stage3.df), ]
rs.seedling_for.late.veget.all <- rownames_to_column(rs.seedling_for.late.veget.all,var = "asvID")
dim(rs.seedling_for.late.veget.all) # 16
# - Early vegetative
rs.early.veget_for.late.veget.all <- rsphere_stage2_nm.rel.indic.df[rownames(rsphere_stage2_nm.rel.indic.df) %in%
                                                                      rownames(rsphere_INdSN.sharedASVs_stage3.df), ]
rs.early.veget_for.late.veget.all <- rownames_to_column(rs.early.veget_for.late.veget.all,var = "asvID")
dim(rs.early.veget_for.late.veget.all) # 13
# - Early reproductive
rs.early.reprod_for.late.veget.all <- rsphere_stage4_nm.rel.indic.df[rownames(rsphere_stage4_nm.rel.indic.df) %in%
                                                                       rownames(rsphere_INdSN.sharedASVs_stage3.df), ]
rs.early.reprod_for.late.veget.all <- rownames_to_column(rs.early.reprod_for.late.veget.all,var = "asvID")
dim(rs.early.reprod_for.late.veget.all) # 16
# - Late reproductive
rs.late.reprod_for.late.veget.all <- rsphere_stage5_nm.rel.indic.df[rownames(rsphere_stage5_nm.rel.indic.df) %in%
                                                                      rownames(rsphere_INdSN.sharedASVs_stage3.df), ]
rs.late.reprod_for.late.veget.all <- rownames_to_column(rs.late.reprod_for.late.veget.all,var = "asvID")
dim(rs.late.reprod_for.late.veget.all) # 16
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rsphere_INdSN.sharedASVs_stage3.df <- rownames_to_column(rsphere_INdSN.sharedASVs_stage3.df, var = "asvID")
rs.late.veget.selected.all <- bind_rows(rs.seedling_for.late.veget.all, rs.early.veget_for.late.veget.all, rsphere_INdSN.sharedASVs_stage3.df,
                                        rs.early.reprod_for.late.veget.all, rs.late.reprod_for.late.veget.all)
View(rs.late.veget.selected.all)

# 4.) Early Reproductive
# Load the indicator taxa which fall in the above prediction in the Neutral model
rsphere_INdSN.sharedASVs_stage4 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizosphere_selected_indicator_early.reproductive.csv")
dim(rsphere_INdSN.sharedASVs_stage4) # 14
# Re-name the old taxonomy SILVA 138 column name
rsphere_INdSN.sharedASVs_stage4 <- rsphere_INdSN.sharedASVs_stage4 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_INdSN.sharedASVs_stage4 <- rsphere_INdSN.sharedASVs_stage4 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_INdSN.sharedASVs_stage4 <- rsphere_INdSN.sharedASVs_stage4 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_INdSN.sharedASVs_stage4.df <- rsphere_INdSN.sharedASVs_stage4 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_INdSN.sharedASVs_stage4.df$point_class <- ifelse(rsphere_INdSN.sharedASVs_stage4.df$above_CI == TRUE,
                                                         "Above",ifelse(rsphere_INdSN.sharedASVs_stage4.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_INdSN.sharedASVs_stage4.df$point_class <- factor(rsphere_INdSN.sharedASVs_stage4.df$point_class,
                                                         levels = c("Above", "Neutral", "Below"))
rsphere_INdSN.sharedASVs_stage4.df$Stage <- "Early Reproductive"
rsphere_INdSN.sharedASVs_stage4.df$stage_aggregate <- as.factor(rsphere_INdSN.sharedASVs_stage4.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_INdSN.sharedASVs_stage4.df <- rsphere_INdSN.sharedASVs_stage4.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rsphere_INdSN.sharedASVs_stage4.df <- column_to_rownames(rsphere_INdSN.sharedASVs_stage4.df, var = "asvID")
# Subset those Early Reproductive selected ASVs from the other stages:
# - Seedling
rs.seedling_for.early.reprod.all <- rsphere_stage1_nm.rel.indic.df[rownames(rsphere_stage1_nm.rel.indic.df) %in%
                                                                     rownames(rsphere_INdSN.sharedASVs_stage4.df), ]
rs.seedling_for.early.reprod.all <- rownames_to_column(rs.seedling_for.early.reprod.all,var = "asvID")
dim(rs.seedling_for.early.reprod.all) # 10
# - Early vegetative
rs.early.veget_for.early.reprod.all <- rsphere_stage2_nm.rel.indic.df[rownames(rsphere_stage2_nm.rel.indic.df) %in%
                                                                        rownames(rsphere_INdSN.sharedASVs_stage4.df), ]
rs.early.veget_for.early.reprod.all <- rownames_to_column(rs.early.veget_for.early.reprod.all,var = "asvID")
dim(rs.early.veget_for.early.reprod.all) # 13
# - Late vegetative
rs.late.veget_for.early.reprod.all <- rsphere_stage3_nm.rel.indic.df[rownames(rsphere_stage3_nm.rel.indic.df) %in%
                                                                       rownames(rsphere_INdSN.sharedASVs_stage4.df), ]
rs.late.veget_for.early.reprod.all <- rownames_to_column(rs.late.veget_for.early.reprod.all,var = "asvID")
dim(rs.late.veget_for.early.reprod.all) # 13
# - Late reproductive
rs.late.reprod_for.early.reprod.all <- rsphere_stage5_nm.rel.indic.df[rownames(rsphere_stage5_nm.rel.indic.df) %in%
                                                                        rownames(rsphere_INdSN.sharedASVs_stage4.df), ]
rs.late.reprod_for.early.reprod.all <- rownames_to_column(rs.late.reprod_for.early.reprod.all,var = "asvID")
dim(rs.late.reprod_for.early.reprod.all) # 11
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rsphere_INdSN.sharedASVs_stage4.df <- rownames_to_column(rsphere_INdSN.sharedASVs_stage4.df, var = "asvID")
rs.early.reprod.selected.all <- bind_rows(rs.seedling_for.early.reprod.all, rs.early.veget_for.early.reprod.all, rs.late.veget_for.early.reprod.all,
                                          rsphere_INdSN.sharedASVs_stage4.df, rs.late.reprod_for.early.reprod.all)
View(rs.early.reprod.selected.all)

# 5.) Late Reproductive
# Load the indicator taxa which fall in the above prediction in the Neutral model
rsphere_INdSN.sharedASVs_stage5 <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_Data/Sloan_Indicator_selected/rhizosphere_selected_indicator_late.reproductive.csv")
dim(rsphere_INdSN.sharedASVs_stage5) # 38
# Re-name the old taxonomy SILVA 138 column name
rsphere_INdSN.sharedASVs_stage5 <- rsphere_INdSN.sharedASVs_stage5 %>%
  dplyr::rename(Old.Domain = predictions.Taxon,
                Old.Phylum = predictions.Phylum,
                Old.Class = predictions.Class,
                Old.Order = predictions.Order,
                Old.Family = predictions.Family,
                Old.Genus = predictions.Genus,
                Old.Species = predictions.Species)
# Add Uniform and Weighted Taxonomy from SILVA 144 in the data frame
rsphere_INdSN.sharedASVs_stage5 <- rsphere_INdSN.sharedASVs_stage5 %>%
  left_join(new.tax.unif %>%
              select(asvID, Unif.Domain, Unif.Kingdom, Unif.Phylum, Unif.Class,
                     Unif.Order, Unif.Family, Unif.Genus), by = "asvID") %>%
  left_join(new.tax.weight %>%
              select(asvID, Weight.Domain, Weight.Kingdom, Weight.Phylum, Weight.Class,
                     Weight.Order, Weight.Family, Weight.Genus, Weight.Species), by = "asvID")
# Relocate the new taxonomy columns
rsphere_INdSN.sharedASVs_stage5 <- rsphere_INdSN.sharedASVs_stage5 %>%
  relocate(Unif.Domain, Weight.Domain, .after = Old.Domain) %>%
  relocate(Unif.Kingdom, Weight.Kingdom, .after = Weight.Domain) %>%
  relocate(Unif.Phylum, Weight.Phylum, .after = Old.Phylum) %>%
  relocate(Unif.Class, Weight.Class, .after = Old.Class) %>%
  relocate(Unif.Order, Weight.Order, .after = Old.Order) %>%
  relocate(Unif.Family, Weight.Family, .after = Old.Family) %>%
  relocate(Unif.Genus, Weight.Genus, .after = Old.Genus) %>%
  relocate(Weight.Species, .after = Old.Species)
# Make some edits
rsphere_INdSN.sharedASVs_stage5.df <- rsphere_INdSN.sharedASVs_stage5 %>%
  mutate(Unif.Phylum = str_remove(Unif.Phylum, "^p__")) %>%
  mutate(Weight.Phylum = str_remove(Weight.Phylum, "^p__"))
rsphere_INdSN.sharedASVs_stage5.df$point_class <- ifelse(rsphere_INdSN.sharedASVs_stage5.df$above_CI == TRUE,
                                                         "Above",ifelse(rsphere_INdSN.sharedASVs_stage5.df$predictions.fit_class == "Below prediction","Below","Neutral"))
rsphere_INdSN.sharedASVs_stage5.df$point_class <- factor(rsphere_INdSN.sharedASVs_stage5.df$point_class,
                                                         levels = c("Above", "Neutral", "Below"))
rsphere_INdSN.sharedASVs_stage5.df$Stage <- "Late Reproductive"
rsphere_INdSN.sharedASVs_stage5.df$stage_aggregate <- as.factor(rsphere_INdSN.sharedASVs_stage5.df$stage_aggregate)
# Calculate mean relative abundance in percent 
rsphere_INdSN.sharedASVs_stage5.df <- rsphere_INdSN.sharedASVs_stage5.df %>%
  dplyr::mutate(mean_relabund_percent = mean_relabund * 100) %>%
  dplyr::relocate(mean_relabund_percent, .after = mean_relabund)
rsphere_INdSN.sharedASVs_stage5.df <- column_to_rownames(rsphere_INdSN.sharedASVs_stage5.df, var = "asvID")
# Subset those Late Reproductive selected ASVs from the other stages:
# - Seedling
rs.seedling_for.late.reprod.all <- rsphere_stage1_nm.rel.indic.df[rownames(rsphere_stage1_nm.rel.indic.df) %in%
                                                                    rownames(rsphere_INdSN.sharedASVs_stage5.df), ]
rs.seedling_for.late.reprod.all <- rownames_to_column(rs.seedling_for.late.reprod.all,var = "asvID")
dim(rs.seedling_for.late.reprod.all) # 36
# - Early vegetative
rs.early.veget_for.late.reprod.all <- rsphere_stage2_nm.rel.indic.df[rownames(rsphere_stage2_nm.rel.indic.df) %in%
                                                                       rownames(rsphere_INdSN.sharedASVs_stage5.df), ]
rs.early.veget_for.late.reprod.all <- rownames_to_column(rs.early.veget_for.late.reprod.all,var = "asvID")
dim(rs.early.veget_for.late.reprod.all) # 34
# - Late vegetative
rs.late.veget_for.late.reprod.all <- rsphere_stage3_nm.rel.indic.df[rownames(rsphere_stage3_nm.rel.indic.df) %in%
                                                                      rownames(rsphere_INdSN.sharedASVs_stage5.df), ]
rs.late.veget_for.late.reprod.all <- rownames_to_column(rs.late.veget_for.late.reprod.all,var = "asvID")
dim(rs.late.veget_for.late.reprod.all) # 34
# - Early reproductive
rs.early.reprod_for.late.reprod.all <- rsphere_stage4_nm.rel.indic.df[rownames(rsphere_stage4_nm.rel.indic.df) %in%
                                                                        rownames(rsphere_INdSN.sharedASVs_stage5.df), ]
rs.early.reprod_for.late.reprod.all <- rownames_to_column(rs.early.reprod_for.late.reprod.all,var = "asvID")
dim(rs.early.reprod_for.late.reprod.all) # 37
# Combine those selected Early Vegetative ASVs with the same ASVs from other stages
rsphere_INdSN.sharedASVs_stage5.df <- rownames_to_column(rsphere_INdSN.sharedASVs_stage5.df, var = "asvID")
rs.late.reprod.selected.all <- bind_rows(rs.seedling_for.late.reprod.all, rs.early.veget_for.late.reprod.all, rs.late.veget_for.late.reprod.all,
                                         rs.early.reprod_for.late.reprod.all, rsphere_INdSN.sharedASVs_stage5.df)
View(rs.late.reprod.selected.all)

# Combine all selected indicator ASVs of certain stage with the same ASVs from other stages
rs.indic.above.all <- bind_rows(rs.seedling.selected.all, rs.early.veget.selected.all, rs.late.veget.selected.all,
                                rs.early.reprod.selected.all, rs.late.reprod.selected.all)
# Save the data
#write.csv(rs.indic.above.all, file="/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/SloanNeutral_Fina/Rdata/Chamber_Study/Tidy_R_Figures/Fina_R_data/rhizosphere_indic_above_all.csv")
View(rs.indic.above.all)




































