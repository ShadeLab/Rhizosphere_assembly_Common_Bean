title: "Rhizo Assembly Plant Health Analysis"
author: "Ari Fina Bintarti"
date: "2024-02-01"

library(ggplot2)
library(ggtext)
library(rstatix)

# Summary
# Plants were grown with or without zebularine treatment, labeled as "Zeb" and "Control".
# The plants treated with zeb grew visibly slower than the control plants. I harvested plant health measurements of plant height and above ground biomass to track growth progression. Below ground biomass was not collected, as the full root system of each plant was frozen for microbiome analysis. Photos of each root system were taken for analysis with imageJ. 
# In the height and biomass data, I will be assessing if the plants harvested by growth stage in zeb and control were the same size, even though they were harvested on different days. I will also be assessing if the plants harvested by time are statistically different due to the different growth rates between zeb and control.

# Load the data
setwd("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/PlantBiomass_Fina/")
data_rhizo_assembly <- read.csv("rhizosphere_assembly_plant_data_ed.csv", header=TRUE)
str(data_rhizo_assembly)

# Change character and integer into factor
data_rhizo_assembly$nodule <- as.numeric(data_rhizo_assembly$nodule)
data_rhizo_assembly[sapply(data_rhizo_assembly, is.integer)] <- lapply(data_rhizo_assembly[sapply(data_rhizo_assembly, is.integer)], as.factor)
data_rhizo_assembly[sapply(data_rhizo_assembly, is.character)] <- lapply(data_rhizo_assembly[sapply(data_rhizo_assembly, is.character)], as.factor)
str(data_rhizo_assembly)

# Subset into control vs. zeb growth and control vs. zeb time, then by time point
time_data <- subset(data_rhizo_assembly, series == "time")
time_3 <- subset(time_data, stage == "day3")
time_7 <- subset(time_data, stage == "day7")
time_14 <- subset(time_data, stage == "day14")
time_21 <- subset(time_data, stage == "day21")
time_35 <- subset(time_data, stage == "day35")
time_49 <- subset(time_data, stage == "day49")
time_63 <- subset(time_data, stage == "day63")

growth_data <- subset(data_rhizo_assembly, series == "growth")
growth_V1 <- subset(growth_data, stage == "V1")
growth_V2 <- subset(growth_data, stage == "V2")
growth_V3 <- subset(growth_data, stage == "V3")
growth_V4 <- subset(growth_data, stage == "V4")
growth_R1 <- subset(growth_data, stage == "R1")
growth_R4 <- subset(growth_data, stage == "R4")
growth_R7 <- subset(growth_data, stage == "R7")

### Statistics

# 1. Dry weight: Time Series

# 
# Compare dry weight between control vs zeb at each time point using Welch two sample t-test 
dry_time_7 <- t.test(postdry ~ group, data = time_7, na.action = "na.omit")
dry_time_7 # ns p-value = 0.05
dry_time_14 <- t.test(postdry ~ group, data = time_14, na.action = "na.omit")
dry_time_14 # * p-value = 0.03
dry_time_21 <- t.test(postdry ~ group, data = time_21, na.action = "na.omit")
dry_time_21 # *** p-value = 0.0009
dry_time_35 <- t.test(postdry ~ group, data = time_35, na.action = "na.omit")
dry_time_35 # **** p-value = 6.352e-05
dry_time_49 <- t.test(postdry ~ group, data = time_49, na.action = "na.omit")
dry_time_49 # **** p-value = 2.826e-05
dry_time_63 <- t.test(postdry ~ group, data = time_63, na.action = "na.omit")
dry_time_63 # *** p-value = 0.0001

# 2. Dry weight: Growth Series

# Compare dry weight between control vs zeb at each growth stage using Welch two sample t-test
dry_growth_V1 <- t.test(postdry ~ group, data = growth_V1, na.action = "na.omit")
dry_growth_V1 # ns p-value = 0.288
dry_growth_V2 <- t.test(postdry ~ group, data = growth_V2, na.action = "na.omit")
dry_growth_V2 # ** p-value = 0.004
dry_growth_V3 <- t.test(postdry ~ group, data = growth_V3, na.action = "na.omit")
dry_growth_V3 # ** p-value = 0.005
dry_growth_V4 <- t.test(postdry ~ group, data = growth_V4, na.action = "na.omit")
dry_growth_V4 # *** p-value = 0.0001
dry_growth_R1 <- t.test(postdry ~ group, data = growth_R1, na.action = "na.omit")
dry_growth_R1 # * p-value = 0.01
dry_growth_R4 <- t.test(postdry ~ group, data = growth_R4, na.action = "na.omit")
dry_growth_R4 # ns p-value = 0.37
dry_growth_R7 <- t.test(postdry ~ group, data = growth_R7, na.action = "na.omit")
dry_growth_R7 # ns p-value = 0.25

# 3. Root length: Time Series

# Compare root length between control vs zeb at each time point using Welch two sample t-test
root_length_3 <- t.test(mean_root_length_mm ~ group, data = time_3, na.action = "na.omit")
root_length_3 # ** p-value = 0.003
root_length_7 <- t.test(mean_root_length_mm ~ group, data = time_7, na.action = "na.omit")
root_length_7 # **** p-value = 1.729e-05
root_length_14 <- t.test(mean_root_length_mm ~ group, data = time_14, na.action = "na.omit")
root_length_14 # **** p-value = 1.086e-05
root_length_21 <- t.test(mean_root_length_mm ~ group, data = time_21, na.action = "na.omit")
root_length_21 # **** p-value = 2.992e-06
root_length_35 <- t.test(mean_root_length_mm ~ group, data = time_35, na.action = "na.omit")
root_length_35 # ** p-value = 0.001
root_length_49 <- t.test(mean_root_length_mm ~ group, data = time_49, na.action = "na.omit")
root_length_49 # ns p-value = 0.06
root_length_63 <- t.test(mean_root_length_mm ~ group, data = time_63, na.action = "na.omit")
root_length_63 # ** p-value = 0.001

# 4. Root length: Growth Series

# Compare root length between control vs zeb at each growth stage using Welch two sample t-test
root_length_V1 <- t.test(mean_root_length_mm ~ group, data = growth_V1, na.action = "na.omit")
root_length_V1 # **** p-value = 5.836e-05
root_length_V2 <- t.test(mean_root_length_mm ~ group, data = growth_V2, na.action = "na.omit")
root_length_V2 # **** p-value = 2.464e-05
root_length_V3 <- t.test(mean_root_length_mm ~ group, data = growth_V3, na.action = "na.omit")
root_length_V3 # *** p-value = 0.0001
root_length_V4 <- t.test(mean_root_length_mm ~ group, data = growth_V4, na.action = "na.omit")
root_length_V4 # ns p-value = 0.5
root_length_R1 <- t.test(mean_root_length_mm ~ group, data = growth_R1, na.action = "na.omit")
root_length_R1 # ns p-value = 0.06
root_length_R4 <- t.test(mean_root_length_mm ~ group, data = growth_R4, na.action = "na.omit")
root_length_R4 # * p-value = 0.01
root_length_R7 <- t.test(mean_root_length_mm ~ group, data = growth_R7, na.action = "na.omit")
root_length_R7 # * p-value = 0.04

# 5. Nodule count: Time Series

# Compare nodule count between control vs zeb at each time point using Welch two sample t-test
nodules_3 <- t.test(nodule ~ group, data = time_3, na.action = "na.omit")
nodules_3 # NA
nodules_7 <- t.test(nodule ~ group, data = time_7, na.action = "na.omit")
nodules_7 # ns
nodules_14 <- t.test(nodule ~ group, data = time_14, na.action = "na.omit")
nodules_14 # **
nodules_21 <- t.test(nodule ~ group, data = time_21, na.action = "na.omit")
nodules_21 # **
nodules_35 <- t.test(nodule ~ group, data = time_35, na.action = "na.omit")
nodules_35 # ****
nodules_49 <- t.test(nodule ~ group, data = time_49, na.action = "na.omit")
nodules_49 # ns
nodules_63 <- t.test(nodule ~ group, data = time_63, na.action = "na.omit")
nodules_63 # ns

# 6. Nodule count: Growth Series

# Compare nodule count between control vs zeb at each growth stage using Welch two sample t-test
nodules_V1 <- t.test(nodule ~ group, data = growth_V1, na.action = "na.omit")
nodules_V1 # ns .36
nodules_V2 <- t.test(nodule ~ group, data = growth_V2, na.action = "na.omit")
nodules_V2 # ns .54
nodules_V3 <- t.test(nodule ~ group, data = growth_V3, na.action = "na.omit")
nodules_V3 # **
nodules_V4 <- t.test(nodule ~ group, data = growth_V4, na.action = "na.omit")
nodules_V4 # **
nodules_R1 <- t.test(nodule ~ group, data = growth_R1, na.action = "na.omit")
nodules_R1 # *
nodules_R4 <- t.test(nodule ~ group, data = growth_R4, na.action = "na.omit")
nodules_R4 # ns .68
nodules_R7 <- t.test(nodule ~ group, data = growth_R7, na.action = "na.omit")
nodules_R7 # ns .26

# tidy up and summarize the data
data_rhizo_assembly.tdy <- data_rhizo_assembly %>%
  group_by(group, series,group.series,stage,stage_all_ed,age_days,PlantID.2) %>%
  summarize(N = n(),Mean.height=mean(height),
            Mean.postdry=if(all(is.na(postdry))) NA_real_ else mean(postdry, na.rm = T),
            Mean.root=if(all(is.na(mean_root_length_mm))) NA_real_ else mean(mean_root_length_mm, na.rm = T),
            Mean.nod=if(all(is.na(nodule))) NA_real_ else mean(nodule, na.rm = T),
            Mean.cp.rhizos=if(all(is.na(rhizos_cp_num))) NA_real_ else mean(rhizos_cp_num, na.rm = T),
            Mean.cp.rhizop=if(all(is.na(rhizop_cp_num))) NA_real_ else mean(rhizop_cp_num, na.rm = T),
            SD.height=sd(height),
            SD.postdry=if(all(is.na(postdry))) NA_real_ else sd(postdry, na.rm = T),
            SD.root=if(all(is.na(mean_root_length_mm))) NA_real_ else sd(mean_root_length_mm, na.rm = T),
            SD.nod=if(all(is.na(nodule))) NA_real_ else sd(nodule, na.rm = T),
            SD.cp.rhizos=if(all(is.na(rhizos_cp_num))) NA_real_ else sd(rhizos_cp_num, na.rm = T),
            SD.cp.rhizop=if(all(is.na(rhizop_cp_num))) NA_real_ else sd(rhizop_cp_num, na.rm = T),
            ymin.height=Mean.height-SD.height,
            ymax.height=Mean.height+SD.height,
            ymin.postdry=Mean.postdry-SD.postdry,
            ymax.postdry=Mean.postdry+SD.postdry,
            SE.postdry.low = Mean.postdry - (SD.postdry/sqrt(N)),
            SE.postdry.high = Mean.postdry + (SD.postdry/sqrt(N)),
            SE.root.low = Mean.root - (SD.root/sqrt(N)),
            SE.root.high = Mean.root + (SD.root/sqrt(N)),
            SE.nod.low = Mean.nod - (SD.nod/sqrt(N)),
            SE.nod.high = Mean.nod + (SD.nod/sqrt(N)),
            SE.cp.rhizos.low = Mean.cp.rhizos - (SD.cp.rhizos/sqrt(N)),
            SE.cp.rhizos.high = Mean.cp.rhizos + (SD.cp.rhizos/sqrt(N)),
            SE.cp.rhizop.low = Mean.cp.rhizop - (SD.cp.rhizop/sqrt(N)),
            SE.cp.rhizop.high = Mean.cp.rhizop + (SD.cp.rhizop/sqrt(N)))
View(data_rhizo_assembly.tdy)

# OR directly read the saved dataframe below:
data_rhizo_assembly.tdy.df <- read.csv("/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/PlantBiomass_Fina/data_rhizo_assembly.tdy.df.csv", sep=",", row.names = 1)
View(data_rhizo_assembly.tdy.df)
# rename
stage_all_order <- c("PE", "VE","VC","V1","V2","V3", "V4", 
                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "RH")
data_rhizo_assembly.tdy.df$stage_all_ed <- factor(data_rhizo_assembly.tdy.df$stage_all_ed, levels = c("Pre_emergence", "VE","VC",
                                                                                        "V1","V2","V3","V4","R1","R2",
                                                                                        "R3","R4","R5","R6","R7","RH"),
                                           labels = c("PE", "VE","VC",
                                                      "V1","V2","V3","V4","R1","R2",
                                                      "R3","R4","R5","R6","R7","RH"))


# Making Plots

# define the x-axis time scales for the control and delay treatment
ctr.t <- c(" ","3","7"," ","14", "21"," "," "," ", "35"," "," ", "49"," ", "63")
zeb.t <- c("3","7","14","21"," ", "35"," "," ","49"," "," ", "63"," ", " "," ")

# assign the x-axis color to match the legend
library(grid)
ctr.day <- textGrob("Control", gp=gpar(fontsize=16,col="#0072B2"))
ctr.3 <- textGrob("3", gp=gpar(fontsize=16,col="#0072B2"))
ctr.7 <- textGrob("7", gp=gpar(fontsize=16,col="#0072B2"))
ctr.14 <- textGrob("14", gp=gpar(fontsize=16,col="#0072B2"))
ctr.21 <- textGrob("21", gp=gpar(fontsize=16,col="#0072B2"))
ctr.35 <- textGrob("35", gp=gpar(fontsize=16,col="#0072B2"))
ctr.49 <- textGrob("49", gp=gpar(fontsize=16,col="#0072B2"))
ctr.63 <- textGrob("63", gp=gpar(fontsize=16,col="#0072B2"))

zeb.day <- textGrob("Delayed", gp=gpar(fontsize=16,col="#D55E00"))
zeb.3 <- textGrob("3", gp=gpar(fontsize=16,col="#D55E00"))
zeb.7 <- textGrob("7", gp=gpar(fontsize=16,col="#D55E00"))
zeb.14 <- textGrob("14", gp=gpar(fontsize=16,col="#D55E00"))
zeb.21 <- textGrob("21", gp=gpar(fontsize=16,col="#D55E00"))
zeb.35 <- textGrob("35", gp=gpar(fontsize=16,col="#D55E00"))
zeb.49 <- textGrob("49", gp=gpar(fontsize=16,col="#D55E00"))
zeb.63 <- textGrob("63", gp=gpar(fontsize=16,col="#D55E00"))
# make border for the days in the bottom of panel
border <- rectGrob(x = unit(0.44, "npc"), #0.46
                   y = unit(-0.28, "npc"), #-0.16
                   width = unit(1.11, "npc"), #1.08
                   height = unit(0.15, "npc"),
                   #just = "centre",
                   gp = gpar(lty=2,lwd = 1, col = "black", fill="#00000000"))
# border title
anno.title <- textGrob("Time\n(days)", gp=gpar(fontsize=14), rot=90)
# assign colors
group.series.col <- c("#0072B2","#D55E00")

# PLANT ABOVE DRY BIOMASS

# Dry weight stat annotation for growth series
dry_anno_LC_growth <- data.frame(xstar = c(3.72, 4.85, 5.85, 6.72, 8, 11.2, 14.2), 
                                 ystar = c(-0.25, -0.45, -0.45, -0.45, -0.45, -0.25, -0.25), 
                                 lab = c( "ns", "**", "**", "***", "*", "ns", "ns"))

# Dry weight stat annotation for time series
dry_anno_LC_time <- data.frame(x1 = c(2, 3, 4,6,9,12), x2 = c(3, 5, 6,10,13,15), 
                               y1 = c(1, 1.4, 2, 5.3, 10.6,12.5), y2 = c(1.2, 1.6, 2.2, 5.5, 10.8, 12.7), 
                               xstar = c(2.5, 4, 5, 8, 11, 13.5), ystar = c(1.6, 1.7, 2.3, 5.6, 11, 12.8),
                               lab = c( "ns", "*", "***", "****", "****", "***"))
# Rename
data_rhizo_assembly.tdy.df$growth.series <- "Growth Series"

# Plot: Above-ground dry biomass - Growth Series
PB_LC_growth <- ggplot(data_rhizo_assembly.tdy.df, 
                       aes(x=factor(stage_all_ed, level=stage_all_order), 
                           y=Mean.postdry, colour = group, shape = series)) +
  geom_point(aes(group=group.series,colour = group), size=3, alpha=0.4)+ #data=data_rhizo_assembly.tdy.df[data_rhizo_assembly.tdy.df$series=="growth",])+ 
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.postdry.low), ymax=ifelse(series == "time", NA, SE.postdry.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ #data=data_rhizo_assembly.tdy.df[data_rhizo_assembly.tdy.df$series=="growth",])+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ #data=data_rhizo_assembly.tdy.df[data_rhizo_assembly.tdy.df$series=="growth",]) +
  theme_bw() +
  facet_grid(. ~ growth.series)+
  scale_y_continuous(limits = c(-0.5, 13))+
  ylab("Above-ground Plant Biomass (g)")+
  labs(title="A", x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "solid", "time" = "blank")) +
  scale_shape_manual(values =c("growth" = 16, "time" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position = c(x=0.1, y=0.8),
        legend.box.background = element_rect(color = "black", linewidth = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.6,xmax=-0.2,ymin=-7.8,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-7.8,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-7.8,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-7.8,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-7.8,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-7.8,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-7.8,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-7.8,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-1.3,xmax=-0.2,ymin=-9.5,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-9.5,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-9.5,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-9.5,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-9.5,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-9.5,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-9.5,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-9.5,ymax=-0.01)+
  annotation_custom(anno.title,xmin=-3,xmax=-2.5,ymin=-8.5,ymax=-0.01)+ #-1.38
  #annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,1,4), "lines")) + #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = dry_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")

PB_LC_growth


# 2.) Subset the data set Time Series

data_rhizo_assembly.tdy.df$time.series <- "Time Series"
PB_LC_time <- ggplot(data_rhizo_assembly.tdy.df, 
                     aes(x=factor(stage_all_ed, level=stage_all_order), 
                         y=Mean.postdry, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.postdry.low), ymax=ifelse(series == "growth", NA, SE.postdry.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(-0.5, 13))+
  facet_grid(. ~ time.series)+
  ylab("Above-ground Plant Biomass (g)")+
  labs(title="B", x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(labels = c("Pre_emergence" = "Pre\nemergence"),expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "blank", "time" = "dashed")) +
  scale_shape_manual(values =c("time" = 17, "growth" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.3,xmax=-0.2,ymin=-9.55,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-9.55,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-9.55,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-9.55,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-9.55,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-9.55,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-9.55,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-9.55,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-1.1,xmax=-0.2,ymin=-11.8,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-11.8,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-11.8,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-11.8,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-11.8,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-11.8,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-11.8,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-11.8,ymax=-0.01)+
  annotation_custom(anno.title,xmin=-1.9,xmax=-2.5,ymin=-10.7,ymax=-0.01)+ #-1.38
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,4.5), "lines")) + #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = dry_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  geom_segment(inherit.aes=FALSE, data =dry_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = dry_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = dry_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

PB_LC_time

library(patchwork)
PBLC.combined <- PB_LC_growth / PB_LC_time  #plot_layout(guides = "collect")
PBLC.combined

setwd('/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/PlantBiomass_Fina/')
ggsave("Plant_Above_Dry_Biomass.tiff",
       PBLC.combined , device = "tiff",
       width = 10, height =12, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")

########################################################################################
# ROOT LENGTH

# Mean Root length with standard error

# Mean Root length stat annotation for growth series
root_anno_LC_growth <- data.frame(xstar = c(3.5, 4.6, 5.8, 6.8, 8, 11.05, 14.09), 
                                  ystar = c(0, 0, 0, 0, 0, 0, 0), 
                                  lab = c( "****", "****", "***", "ns", "ns", "*", "*"))
root_anno_LC_growth
# Mean Root length stat annotation for time series
root_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4, 6, 9, 12), x2 = c(2, 3, 5, 6, 10, 13, 15), 
                                y1 = c(80, 165, 350, 400, 410,420, 400), y2 = c(82, 167, 352, 402, 412, 422, 402), 
                                xstar = c(1.5, 2.5, 4, 5, 8, 11,13.5), ystar = c(83, 169, 354, 403, 413, 440,403),
                                lab = c( "**", "****", "****", "****", "**", "ns","**"))
root_anno_LC_time


# 1. Root length - Growth Series

data_rhizo_assembly.tdy.df$rl.growth.series <- "A1. Growth Series"

RL_LC_growth <- ggplot(data_rhizo_assembly.tdy.df, 
                       aes(x=factor(stage_all_ed, level=stage_all_order), 
                           y=Mean.root, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.root.low), ymax=ifelse(series == "time", NA, SE.root.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ rl.growth.series)+
  scale_y_continuous(limits = c(0, 450))+
  ylab("Root Length (mm)")+
  labs(title="A. Root Length", x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "solid", "time" = "blank")) +
  scale_shape_manual(values =c("growth" = 16, "time" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position = c(x=0.13, y=0.8),
        legend.box.background = element_rect(color = "black", linewidth = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))+  #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = root_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=6,vjust="inward",hjust="inward")


RL_LC_growth

# 2. Root length - Time Series

data_rhizo_assembly.tdy.df$rl.time.series <- "A2. Time Series"

RL_LC_time <- ggplot(data_rhizo_assembly.tdy.df, 
                     aes(x=factor(stage_all_ed, level=stage_all_order), 
                         y=Mean.root, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.root.low), ymax=ifelse(series == "growth", NA, SE.root.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(0, 450))+
  facet_grid(. ~ rl.time.series)+
  ylab("Root Length (mm)")+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(labels = c("Pre_emergence" = "Pre\nemergence"),
                   expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "blank", "time" = "dashed")) +
  scale_shape_manual(values =c("time" = 17, "growth" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.2,ymin=-220,ymax=-70) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-220,ymax=-70) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-220,ymax=-70)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-220,ymax=-70)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-220,ymax=-70)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-220,ymax=-70)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-220,ymax=-70)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-220,ymax=-70)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=-270,ymax=-90) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-270,ymax=-90) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-270,ymax=-90)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-270,ymax=-90)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-270,ymax=-90)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-270,ymax=-90)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-270,ymax=-90)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-270,ymax=-90)+
  annotation_custom(anno.title,xmin=-2.4,xmax=-2,ymin=-250,ymax=-80)+
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) + #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = root_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  geom_segment(inherit.aes=FALSE, data =root_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = root_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = root_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

RL_LC_time


RLLC.combined <- RL_LC_growth / RL_LC_time  #plot_layout(guides = "collect")
RLLC.combined

########################################################################################
# NODULE COUNT

# Nodule count stat annotation for growth series
nod_anno_LC_growth <- data.frame(xstar = c(3.8, 4.8, 5.9, 6.9, 8, 11.15, 14.13), 
                                 ystar = c(-2, -2, -2, -2, -2, -2, -2), 
                                 lab = c( "ns", "ns", "**", "**", "*", "ns", "ns"))
nod_anno_LC_growth
# Nodule count stat annotation for time series
nod_anno_LC_time <- data.frame(x1 = c(2, 3, 4, 6, 9, 12), x2 = c(3, 5, 6, 10, 13, 15), 
                               y1 = c(10, 50, 148, 200,210, 125), y2 = c(12, 52, 150, 202, 212, 127), 
                               xstar = c(2.5, 4, 5, 8, 11,13.5), ystar = c(20, 54, 151, 204, 220,132),
                               lab = c("ns", "**", "**", "****", "ns","ns"))
nod_anno_LC_time

# 1. Nodule count - Growth Series

data_rhizo_assembly.tdy.df$nod.growth.series <- "B1. Growth Series"

NOD_LC_growth <- ggplot(data_rhizo_assembly.tdy.df, 
                        aes(x=factor(stage_all_ed, level=stage_all_order), 
                            y=Mean.nod, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.nod.low), ymax=ifelse(series == "time", NA, SE.nod.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ nod.growth.series)+
  scale_y_continuous(limits = c(-2, 230))+
  ylab("Nodule Count")+
  labs(title="B. Nodule Count", x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "solid", "time" = "blank")) +
  scale_shape_manual(values =c("growth" = 16, "time" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        legend.box.background = element_rect(color = "black", linewidth = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))+  
  geom_text(inherit.aes=FALSE, data = nod_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=6,vjust="inward",hjust="inward")


NOD_LC_growth

# 2. Nodule count - Time Series

data_rhizo_assembly.tdy.df$nod.time.series <- "B2. Time Series"
NOD_LC_time <- ggplot(data_rhizo_assembly.tdy.df, 
                      aes(x=factor(stage_all_ed, level=stage_all_order), 
                          y=Mean.nod, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.nod.low), ymax=ifelse(series == "growth", NA, SE.nod.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ nod.time.series)+
  scale_y_continuous(limits = c(-2, 230))+
  ylab("Nodule Count")+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(labels = c("Pre_emergence" = "Pre\nemergence"),
                   expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "blank", "time" = "dashed")) +
  scale_shape_manual(values =c("time" = 17, "growth" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.05,ymin=-110,ymax=-43) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-110,ymax=-43) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-110,ymax=-43)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-110,ymax=-43)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-110,ymax=-43)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-110,ymax=-43)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-110,ymax=-43)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-110,ymax=-43)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=-136,ymax=-55) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-136,ymax=-55) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-136,ymax=-55)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-136,ymax=-55)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-136,ymax=-55)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-136,ymax=-55)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-136,ymax=-55)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-136,ymax=-55)+
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) + #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = nod_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  geom_segment(inherit.aes=FALSE, data =nod_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = nod_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = nod_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

NOD_LC_time


NODLC.combined <- NOD_LC_growth / NOD_LC_time  #plot_layout(guides = "collect")
NODLC.combined

# Combine between mean root length and nodule count
RL_NOD.all <- RLLC.combined | NODLC.combined

setwd('/Users/emiliedehon/Nextcloud/Microrescue_share/Projects_Papers/Rhizosphere_assembly_Common_Bean-main/PlantBiomass_Fina/')
ggsave("Root_Length_Nodule_Count.tiff",
       RL_NOD.all , device = "tiff",
       width = 16.5, height =12, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")



