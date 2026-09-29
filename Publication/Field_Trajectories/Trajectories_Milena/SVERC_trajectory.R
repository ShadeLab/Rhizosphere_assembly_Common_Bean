
map=map_final[map_final$ID %in% colnames(otu_final),]
str(map)

map=map[order(map$ID),]

colnames(otu_final)==map$ID

# Core taxa abundance
map <- map[!(map$ID %in% PowersoilMRC),]
otu_final <- otu_final[,!(colnames(otu_final) %in% PowersoilMRC)]
rel.abun <- decostand(otu_final, method = 'total', MARGIN = 2)

global_abun_df<- data.frame(otu=rownames(rel.abun), rel.abun) %>%
  gather(ID, abun, -otu) %>%
  left_join(map) %>%
  filter(otu %in% global_core,
         !is.na(Timepoint)) %>%
  group_by(ID,Timepoint, Site, Compartment) %>%
  summarise(coreNo=sum(abun),
            timeNo=length(unique(ID)),
            rel=coreNo/timeNo) 
str(global_abund_df)

# Load necessary libraries
library(vegan)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Filter data for SVERC
sverc_data <- global_abun_df %>% filter(Site == "SVERC")

# Split data by Compartment
sverc_compartments <- split(sverc_data, sverc_data$Compartment)

# Function to compute ecological trajectory
compute_trajectory <- function(compartment_data) {
  compartment_data <- compartment_data[order(compartment_data$Timepoint), ]
  
  # Compute Bray-Curtis distance
  bray_dist <- vegdist(as.matrix(compartment_data$rel), method = "bray")
  
  # Perform PCoA
  pcoa_result <- cmdscale(bray_dist, k = 2)
  
  # Create results data frame
  pcoa_df <- data.frame(
    ID = compartment_data$ID,
    Timepoint = compartment_data$Timepoint,
    PC1 = pcoa_result[,1],
    PC2 = pcoa_result[,2],
    Compartment = compartment_data$Compartment
  )
  
  # Compute distances
  pcoa_df <- pcoa_df %>%
    arrange(Timepoint) %>%
    mutate(
      dist_to_previous = c(NA, sqrt(diff(PC1)^2 + diff(PC2)^2)),
      cumulative_distance = c(0, cumsum(na.omit(dist_to_previous)))
    )
  
  # Calculate total distance for each trajectory (distance from first to last point)
  total_distance <- sqrt((pcoa_df$PC1[nrow(pcoa_df)] - pcoa_df$PC1[1])^2 + (pcoa_df$PC2[nrow(pcoa_df)] - pcoa_df$PC2[1])^2)
  pcoa_df$total_distance <- total_distance
  
  return(pcoa_df)
}

# Apply function to each compartment
sverc_results <- bind_rows(lapply(sverc_compartments, compute_trajectory))

# Statistical tests
# PERMANOVA (adonis test) for community differences between compartments
perm_result_sverc <- adonis2(global_abun_df$rel ~ global_abun_df$Compartment, data = global_abun_df %>% filter(Site == "SVERC"), method = "bray")

# Wilcoxon test for cumulative distances
wilcox_test_cumulative_sverc <- wilcox.test(cumulative_distance ~ Compartment, data = sverc_results)

# Wilcoxon test for total distances
wilcox_test_total_sverc <- wilcox.test(total_distance ~ Compartment, data = sverc_results)

# Export results to Excel
write.xlsx(list(
  "SVERC Results" = sverc_results,
  "PERMANOVA" = data.frame(perm_result_sverc),
  "Wilcoxon Cumulative" = data.frame(
    statistic = wilcox_test_cumulative_sverc$statistic, 
    p_value = wilcox_test_cumulative_sverc$p.value
  ),
  "Wilcoxon Total" = data.frame(
    statistic = wilcox_test_total_sverc$statistic, 
    p_value = wilcox_test_total_sverc$p.value
  )
), "SVERC_ecological_trajectories.xlsx")

# Plot results
ggplot(sverc_results, aes(x = PC1, y = PC2, color = as.factor(Timepoint), group = Compartment)) +
  geom_point(size = 3) +
  geom_line(aes(group = Compartment), arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
  facet_wrap(~ Compartment) +
  theme_minimal() +
  labs(title = "SVERC Ecological Trajectories", x = "PC1", y = "PC2") +
  theme(legend.position = "bottom")

