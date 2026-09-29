
# Script for Rhizhospere microbial assembly Project
# Xipeng Liu
# 2025-02-05

# COG ####
library(ggplot2)

# read gene table
data <- read.table("data/COG_category_grouped_abundance.txt", header=TRUE, row.names=1, sep="\t")
data <- read.table("data/COG_gene_abundance.txt", header=TRUE, row.names=1, sep="\t")
data <- data[, c(1:27)]

head(data)


## PCOA ####
library(vegan)
library(ggalt)

metadata <- read.table("data/metadata_grouped.txt", header=TRUE, col.names=c("Sample", "Group"))

data_pcoa <- t(data)

all.dist <- vegdist(data_pcoa, method = "bray", binary = F)

all.pcoa <- cmdscale(all.dist, k=3, eig = T)
all_pcoa_points <- as.data.frame(all.pcoa$points)
sum_eig <- sum(all.pcoa$eig)
eig_percent <- round(all.pcoa$eig/sum_eig*100, 1)
colnames(all_pcoa_points) <- paste0("PCoA", 1:3)

# combine pcoa and sample info
all_pcoa_result <- cbind(all_pcoa_points, metadata)

head(all_pcoa_result)

### Adonis statistic
set.seed(1)
all.div <- adonis2(all.dist ~ Group, data = metadata, permutations = 999, method = "bray")
all.div

## legend order
all_pcoa_result$Group <- factor(all_pcoa_result$Group, 
                                levels = c("V2","V5","Flowering",
                                           "Pod_filling","Senescence"))

## plot
p1 <- ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Group, group=Group)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(aes(group=Group), size=5, alpha=0.8) +
  #geom_encircle(aes(group=Group), alpha=0.4, show.legend = F, color="black") +
  labs(title = "COG", size=17)+
  theme_bw()+
  scale_color_manual(values = c("#81af4a","#357933","#9f1214","#25567d","#636363"))+
  scale_fill_manual(values = c("#81af4a","#357933","#9f1214","#25567d","#636363"))+
  scale_shape_manual(values=c(19, 17))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=17), axis.title.y=element_text(size=17),
        axis.text.x =element_text(size=13, color='black'), axis.text.y=element_text(size=13, color='black')) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15))+
  geom_text(x=-0.004, y=0.0045, label="p < 0.001; Adonis", color="black", size=5)

p1

ggsave("PCoA_COG_gene.jpeg", dpi=600, width=18, height=13, units="cm", bg = "white")



## Abundance ####
#library(circlize)
#library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

data <- read.table("COG_category_grouped_abundance.txt", header=TRUE, row.names=1, sep="\t")

abundance_cog <- aggregate(. ~ COG_category_2, data = data, FUN = sum)

rownames(abundance_cog) <- abundance_cog[, 1]

abundance_cog <- abundance_cog[, -1]

# pie plot showing percentage of each category
mean_abundance <- rowMeans(abundance_cog)

mean_abundance_df <- data.frame(
  Category = names(mean_abundance),
  Abundance = mean_abundance
)

#read COG category
cog_tab <- read.csv("COG_category.csv", header=TRUE, row.names=1)

merged_cog <- merge(mean_abundance_df, cog_tab, by.x = "Category", by.y = "Category_1")

merged_cog$Category_with_description <- paste(merged_cog$Category, "-", merged_cog$Group_1)

# top 7 categories
top7_categories <- merged_cog %>%
  arrange(desc(Abundance)) %>%
  head(7) %>%
  pull(Category)

merged_cog$label_display <- ifelse(merged_cog$Category %in% top7_categories, merged_cog$Category, "")

# plot
ggplot(merged_cog, aes(x = "", y = Abundance, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  ggtitle("Gene Category Proportions") +
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 10),
        legend.margin = margin(t = 10), 
        legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +  
  scale_fill_manual(values = rainbow(length(merged_cog$Category)),
                    breaks = merged_cog$Category, 
                    labels = merged_cog$Category_with_description) +
  guides(fill = guide_legend(ncol = 1))+
  geom_text(aes(label = label_display), position = position_stack(vjust = 0.5), color = "black", size = 3)

ggsave("Gene_COG_Category_Proportions.jpeg", dpi=600, width=25, height=14, units="cm", bg = "white")






## Cluster Mfuzz ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mfuzz")

library(Mfuzz)

library(dplyr)
library(ggplot2)

# read gene table with COG annotation; 
data <- read.table("COG_gene_abundance.txt", header=TRUE, row.names=1, sep="\t")

# filtering the unknown functions; 
data <- subset(data, !(COG_category == "S" | COG_category == "-" | COG_category == "R"))



# calculate average for each time point
metadata <- read.table("metadata_grouped.txt", header=TRUE, col.names=c("Sample", "Group"))

group_means <- sapply(unique(metadata$Group), function(group) {
  samples <- metadata$Sample[metadata$Group == group]
  rowMeans(data[, samples, drop=FALSE], na.rm=TRUE)
})

group_means <- group_means[, c("V2", "V5", "Flowering", "Pod_filling", "Senescence")]

#write.table(group_means, file="COG_gene_filtered_abundance_avg.txt", sep="\t", quote=FALSE)


gene_table <- read.table("COG_gene_filtered_abundance_avg.txt", header=TRUE, row.names=1, sep="\t")

# filter gene abundance with 0 in more than 3 stages
gene_table <- gene_table[apply(gene_table, 1, function(row) sum(row == 0) <= 2), ]

# filter abundance with low sd
row_sd <- apply(gene_table, 1, sd)

gene_table <- gene_table[row_sd >= 0.1, ]

write.table(gene_table, file="COG_gene_filtered_abundance_final.txt", sep="\t", quote=FALSE)

gene_table <- read.table("COG_gene_filtered_abundance_final.txt", header=TRUE, row.names=1, sep="\t")


# standardize 
eset <- ExpressionSet(assayData=as.matrix(gene_table))

eset <- standardise(eset)

summary(exprs(eset))

# Estimate for optimal fuzzifier.
m=mestimate(eset)

# Elbow Method to determine number of clusters; only run once
withinerror_values <- sapply(2:12, function(c) {
  set.seed(123)
  mfuzz_result <- mfuzz(eset, c=c, m=m)
  mfuzz_result$withinerror
})

# plot Elbow curve
# If the inflection point is not obvious, it selected based on biological significance.
plot(2:12, withinerror_values, type="b", pch=19, xlab="Number of Clusters", ylab="Withinerror")

write.csv(withinerror_values, "withinerror_values.csv")


# analysis and plot
# We use 8 clusters.
set.seed(123)
mfuzz_result <- mfuzz(eset, c=8, m=m)

# save mfuzz_result
save(mfuzz_result, file = "COG_mfuzz_result.RData")

# load mfuzz_result
load("COG_mfuzz_result.RData")

# plot: yellow or green lines: genes with low values; red and purple: genes with high values.
mfuzz.plot2(eset, cl=mfuzz_result, mfrow=c(2,4),
            xlab="Stage",ylab="Abundance changes",
            centre=TRUE, centre.lwd=1,centre.col="black")


# save results
# extract cluster membership
membership <- mfuzz_result$membership
write.csv(membership, file = "COG_mfuzz_membership.csv")

cluster_labels <- mfuzz_result$cluster
write.csv(cluster_labels, file = "COG_mfuzz_cluster_labels.csv")

cluster_centers <- mfuzz_result$centers
write.csv(cluster_centers, file = "COG_mfuzz_cluster_centers.csv")



# gene in each cluster
size <- mfuzz_result$size
write.csv(cluster_centers, file = "COG_mfuzz_size.csv")

# overlap
overlap <- overlap(mfuzz_result)
Ptmp <- overlap.plot(mfuzz_result, over=overlap, thres=0.1)


# alpha cores with COG category

## gene-COG mapping table
cog_map <- read.table("COG_gene_map_filtered.txt", header=FALSE, sep="\t")

cog_map <- subset(cog_map, !(V2 == "S" | V2 == "-" | V2 == "NA"))

print(colnames(cog_map))

#write.csv(cog_map, "cog_map_filtered.csv")

membership_filtered <- acore(eset, mfuzz_result, min.acore=0.7)

print(colnames(membership_filtered[[1]]))

all_clusters <- list()

for (i in seq_along(membership_filtered)) {
  cluster_data <- membership_filtered[[i]]
  if (nrow(cluster_data) == 0) {
    cat("Skipping empty cluster:", i, "\n")
    next
  }
  cluster_data <- merge(cluster_data, cog_map, by.x = "NAME", by.y = "V1", all.x = TRUE)
  cluster_data$Cluster_ID <- i  
  all_clusters[[i]] <- cluster_data  
}

final_data <- bind_rows(all_clusters)

write.csv(final_data, file = "all_clusters_membership_0.7_COG.csv", row.names = FALSE)
cat("Exported all_clusters_membership_0.7.csv\n")






# Count the functional category proportion of different clusters ####
library(dplyr)
library(ggplot2)
library(readxl) # read excel
library(RColorBrewer)

stage_cluster <- read.csv("all_clusters_membership_0.7_COG.csv", row.names = 1)

stage_cluster <- subset(stage_cluster, !(V2 == "S" | V2 == "-" | V2 == "NA"))


## top 50 mem.ship ####
library(RColorBrewer)

df_top_cluster <- stage_cluster %>%
  group_by(Stage) %>%
  top_n(50, MEM.SHIP) %>%
  ungroup()

# calculate percentage
df_top_cluster <- df_top_cluster %>%
  group_by(Stage, COG_category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = Count / sum(Count) * 100)

## legend order
df_top_cluster$Stage <- factor(df_top_cluster$Stage, 
                                levels = c("V2","V5","Flowering",
                                           "Pod_filling","Senescence"))

df_top_cluster <- df_top_cluster %>% arrange(COG_category)

color_palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(df_top_cluster$COG_category)))

p <- ggplot(df_top_cluster, aes(x = Stage, y = Percentage, fill = COG_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_bw() +
  labs(x = "Stage", y = "Percentage (%)", title = "Top 50 genes' COG category") +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=16), 
        axis.title.y=element_text(size=16),
        axis.text.x =element_text(size=13, color='black'), 
        axis.text.y=element_text(size=12, color='black'),
        strip.text = element_text(size = 10, color = "black")) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = 10), 
        legend.box.margin = margin(10, 10, 10, 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) 

p

ggsave("Fig.2e_Top_50_membership_genes'_COG_category.jpeg", dpi=600, width = 16, height = 14, units = "cm")


## significance test of percentage composition and top 100 heatmap ####
library(ggplot2)
library(vegan)
library(tidyr)

v2_cog_mapping <- stage_cluster %>%
  distinct(V2, COG_category)

# stat
cog_percentage <- stage_cluster %>%
  group_by(Stage, V2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Stage) %>%
  mutate(total_count = sum(count),  
         percentage = count / total_count * 100)

cog_percentage_table <- reshape2::dcast(cog_percentage, Stage ~ V2, value.var = "percentage", fill = 0)

data_long <- pivot_longer(cog_percentage_table, 
                          cols = -Stage, 
                          names_to = "COG", 
                          values_to = "Percentage")

column_means <- colMeans(select(cog_percentage_table, -Stage), na.rm = TRUE)
top_100_categories <- sort(column_means, decreasing = TRUE)[1:100]

top_100_data <- cog_percentage_table %>% 
  select(Stage, any_of(names(top_100_categories)))

data_long_100 <- pivot_longer(top_100_data, 
                              cols = -Stage, 
                              names_to = "COG", 
                              values_to = "Percentage")

# mapping COG_category
data_long_100 <- data_long_100 %>%
  left_join(v2_cog_mapping, by = c("COG" = "V2"))

data_long_100$COG <- reorder(data_long_100$COG, data_long_100$Percentage, FUN = mean, na.rm = TRUE)


## legend order
data_long_100$Stage <- factor(data_long_100$Stage, 
                              levels = c("V2","V5","Flowering",
                                         "Pod_filling","Senescence"))

cog_categories <- unique(data_long_100$COG_category)
num_categories <- length(cog_categories)
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_categories)

ggplot(data_long_100, aes(x = COG, y = Percentage, fill = COG_category)) +
  geom_bar(stat = "identity", size = 2) +
  labs(x = "COG (top 100)", y = "Percentage (%)") +
  labs(title = "Percentage of gene categories", size=17)+
  scale_fill_manual(values = color_palette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(rows = vars(Stage))+
  theme(panel.background=element_rect(fill='white', color='black'), 
        axis.title.x =element_text(size=16), axis.title.y=element_text(size=16),
        axis.text.x =element_text(size=7, color='black'), 
        axis.text.y=element_text(size=8, color='black'),
        strip.text = element_text(size = 10, color = "black")) +
  theme(legend.text=element_text(size=14), legend.title=element_text(size=14))

ggsave("Fig._distribution_gene_categories.jpeg", dpi=600, width = 38, height = 14, units = "cm")


# normal distribution test; p-value < 0.05, suggesting non-normal distribution
by(data_long$Percentage, data_long$Stage, shapiro.test)

# pairwise.wilcox.test
pairwise_result <- pairwise.wilcox.test(data_long$Percentage, data_long$Stage)
print(pairwise_result)


















# KEGG ####
# KEGG shows the similar patthern compared with COG.
library(ggplot2)

# read gene table
data <- read.table("KEGG_ko_grouped_abundance.txt", header=TRUE, row.names=1, sep="\t")

head(data)


## PCOA ####
library(vegan)
library(ggalt)

metadata <- read.table("metadata_grouped.txt", header=TRUE, col.names=c("Sample", "Group"))

data_pcoa <- t(data)

all.dist <- vegdist(data_pcoa, method = "bray", binary = F)

all.pcoa <- cmdscale(all.dist, k=3, eig = T)
all_pcoa_points <- as.data.frame(all.pcoa$points)
sum_eig <- sum(all.pcoa$eig)
eig_percent <- round(all.pcoa$eig/sum_eig*100, 1)
colnames(all_pcoa_points) <- paste0("PCoA", 1:3)

# combine pcoa and sample info
all_pcoa_result <- cbind(all_pcoa_points, metadata)

head(all_pcoa_result)

### Adonis statistic
set.seed(1)
all.div <- adonis2(all.dist ~ Group, data = metadata, permutations = 999, method = "bray")
all.div

## legend order
all_pcoa_result$Group <- factor(all_pcoa_result$Group, 
                                    levels = c("V2","V5","Flowering",
                                               "Pod_filling","Senescence"))

## plot
p1 <- ggplot(all_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Group, group=Group)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep = ""), 
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep = "")) +
  geom_point(aes(group=Group), size=5, alpha=0.8) +
  #geom_encircle(aes(group=Group), alpha=0.4, show.legend = F, color="black") +
  labs(title = "KEGG", size=17)+
  theme_bw()+
  scale_color_manual(values = c("#81af4a","#357933","#9f1214","#25567d","#636363"))+
  scale_fill_manual(values = c("#81af4a","#357933","#9f1214","#25567d","#636363"))+
  scale_shape_manual(values=c(19, 17))+
  theme(panel.background=element_rect(fill='white', color='black'), axis.title.x =element_text(size=17), axis.title.y=element_text(size=17),
        axis.text.x =element_text(size=13, color='black'), axis.text.y=element_text(size=13, color='black')) +
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15))+
  geom_text(x=0.01, y=0.022, label="p < 0.001; Adonis", color="black", size=5)

p1

ggsave("PCoA_KEGG_gene.jpeg", dpi=600, width=18, height=13, units="cm", bg = "white")


