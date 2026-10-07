title: "Rhizo Assembly Plant Health Analysis"
author: "Ari Fina Bintarti"
date: "2024-02-01"

library(ggplot2)
library(ggtext)
library(bestNormalize)


#################################################################################################





#### Load in data
setwd("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Rhizosphere_assembly_Common_Bean/R_Analysis_Files/")
data_rhizo_assembly <- read.csv("rhizosphere_assembly_plant_data_ed.csv", header=TRUE)
str(data_rhizo_assembly)

#change character and integer into factor
data_rhizo_assembly$nodule <- as.numeric(data_rhizo_assembly$nodule)
data_rhizo_assembly[sapply(data_rhizo_assembly, is.integer)] <- lapply(data_rhizo_assembly[sapply(data_rhizo_assembly, is.integer)], as.factor)
data_rhizo_assembly[sapply(data_rhizo_assembly, is.character)] <- lapply(data_rhizo_assembly[sapply(data_rhizo_assembly, is.character)], as.factor)
str(data_rhizo_assembly)
view(data_rhizo_assembly)

# tidy up and summarize the data frame
data_rhizo_assembly$stage_all_ed <- factor(data_rhizo_assembly$stage_all_ed, levels = c("Pre_emergence", "VE","VC",
                                                                                        "V1","V2","V3","V4","R1","R2",
                                                                                        "R3","R4","R5","R6","R7","RH"),
                                           labels = c("PE", "VE","VC",
                                                      "V1","V2","V3","V4","R1","R2",
                                                      "R3","R4","R5","R6","R7","RH"))


# subset data (for stats analysis) into control v zeb growth and control v zeb time, then by time point
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


#____STAT FOR 16S COPIES

# 4. 16S rRNA gene copy number 

# ANOVA - 16S rRNA gene copy number Rhizosphere

# 1.) time series
# between group and age (days) in time series
rhizo.cp.time.aov <- aov(rhizos_cp_num ~ group*age_days, data=time_data)
summary(rhizo.cp.time.aov) # NS
#                Df Sum Sq Mean Sq F value   Pr(>F)
#group           1   17.8  17.840   1.444  0.235
#age_days        6   57.0   9.503   0.769  0.597
#group:age_days  6   80.8  13.472   1.090  0.380

log.rhizo.cp.time.aov <- aov(log.rhizos_cp_num ~ group*age_days, data=time_data)
summary(log.rhizo.cp.time.aov) # NS
#                Df Sum Sq Mean Sq F value   Pr(>F)
#group           1  0.774  0.7745   6.972 0.0107 *
#age_days        6  0.387  0.0645   0.581 0.7441  
#group:age_days  6  1.198  0.1996   1.797 0.1164 

# 2.) growth series
# between group and growth stage in growth series
rhizo.cp.growth.aov <- aov(rhizos_cp_num ~ group*stage, data=growth_data)
summary(rhizo.cp.growth.aov)
#             Df Sum Sq Mean Sq F value   Pr(>F)
#group        1   0.06   0.058   0.016 0.8986  
#stage        6  25.13   4.188   1.179 0.3309  
#group:stage  6  50.49   8.415   2.368 0.0414 *

arcsinh.rhizo.cp.growth.aov <- aov(arcsinh.rhizos_cp_num ~ group*stage, data=growth_data)
summary(arcsinh.rhizo.cp.growth.aov)
#             Df Sum Sq Mean Sq F value   Pr(>F)
#group        1   0.20  0.2039   0.243 0.6243  
#stage        6   8.34  1.3906   1.655 0.1494  
#group:stage  6  13.39  2.2319   2.656 0.0245 *



# check assumption (outliers)
rhizo.cp.time.out <- time_data %>%
  group_by(group) %>%
  identify_outliers(log.rhizos_cp_num) # no outliers
View(rhizo.cp.time.out)

rhizo.cp.growth.out <- growth_data %>%
  group_by(group) %>%
  identify_outliers(rhizos_cp_num) # no outliers
View(rhizo.cp.growth.out)

# Saphiro-Wilk for normality
rhizo.cp.time.SW <- time_data %>%
  group_by(group) %>%
  shapiro_test(log.rhizos_cp_num)
View(rhizo.cp.time.SW) # good

rhizo.cp.growth.SW <- growth_data %>%
  group_by(group) %>%
  shapiro_test(arcsinh.rhizos_cp_num)
View(rhizo.cp.growth.SW) # alright

# Lavene test
rhizo.cp.time.Lave <- time_data %>%
  levene_test(rhizos_cp_num ~ group)
View(rhizo.cp.time.Lave) # good

rhizo.cp.growth.Lave <- growth_data %>%
  levene_test(rhizos_cp_num ~ group)
View(rhizo.cp.growth.Lave) # good


# ANOVA - 16S rRNA gene copy number Rhizoplane

# 1.) time series
# between group and age (days) in time series
rhizop.cp.time.aov <- aov(rhizop_cp_num ~ group*age_days, data=time_data)
summary(rhizop.cp.time.aov) # NS
#                Df Sum Sq Mean Sq F value   Pr(>F)
#ggroup           1   21.5   21.47   1.931 0.170111    
#age_days        6  218.0   36.34   3.269 0.007943 ** 
#group:age_days  6  303.3   50.55   4.547 0.000803 ***

arcsinh.rhizop.cp.time.aov <- aov(arcsinh.rhizop_cp_num ~ group*age_days, data=time_data)
summary(arcsinh.rhizop.cp.time.aov) # NS
#                Df Sum Sq Mean Sq F value   Pr(>F)
#group           1   0.48  0.4835   0.672 0.41580   
#age_days        6  12.41  2.0684   2.875 0.01634 * 
#group:age_days  6  15.82  2.6368   3.665 0.00386 **


# 2.) growth series
# between group and growth stage in growth series
rhizop.cp.growth.aov <- aov(rhizop_cp_num ~ group*stage, data=growth_data)
summary(rhizop.cp.growth.aov) #NS
#             Df Sum Sq Mean Sq F value   Pr(>F)
#group        1   12.9   12.87   0.639  0.427
#stage        6  202.3   33.72   1.675  0.144
#group:stage  6  135.9   22.66   1.125  0.360

arcsinh.rhizop.cp.growth.aov <- aov(arcsinh.rhizop_cp_num ~ group*stage, data=growth_data)
summary(arcsinh.rhizop.cp.growth.aov) #NS
#group        1   0.20  0.1985   0.209  0.650
#stage        6   8.83  1.4712   1.547  0.180
#group:stage  6   6.71  1.1191   1.177  0.332




# check assumption (outliers)
rhizop.cp.time.out <- time_data %>%
  group_by(group) %>%
  identify_outliers(arcsinh.rhizop_cp_num) # no outliers
View(rhizop.cp.time.out)

rhizop.cp.growth.out <- growth_data %>%
  group_by(group) %>%
  identify_outliers(arcsinh.rhizop_cp_num) # no outliers
View(rhizop.cp.growth.out)

# Saphiro-Wilk for normality
rhizop.cp.time.SW <- time_data %>%
  group_by(group) %>%
  shapiro_test(arcsinh.rhizop_cp_num)
View(rhizop.cp.time.SW) # good

rhizop.cp.growth.SW <- growth_data %>%
  group_by(group) %>%
  shapiro_test(arcsinh.rhizop_cp_num)
View(rhizop.cp.growth.SW) # good

# Lavene test
rhizop.cp.time.Lave <- time_data %>%
  levene_test(rhizop_cp_num ~ group)
View(rhizop.cp.time.Lave) # good

rhizop.cp.growth.Lave <- growth_data %>%
  levene_test(rhizop_cp_num ~ group)
View(rhizop.cp.growth.Lave) # good


hist(time_data$rhizos_cp_num)
hist(growth_data$rhizos_cp_num)

hist(time_data$rhizop_cp_num)
hist(growth_data$rhizop_cp_num)

# Transform data

# Rhizosphere
time_data$log.rhizos_cp_num <- log10(time_data$rhizos_cp_num)
hist(time_data$log.rhizos_cp_num)

set.seed(13)
arcsinh.rhizos_cp_num <- arcsinh_x(growth_data$rhizos_cp_num)
growth_data$arcsinh.rhizos_cp_num <- arcsinh.rhizos_cp_num$x.t
hist(growth_data$arcsinh.rhizos_cp_num)

# Rhizoplane
set.seed(13)
arcsinh.rhizop_cp_num <- arcsinh_x(time_data$rhizop_cp_num)
time_data$arcsinh.rhizop_cp_num <- arcsinh.rhizop_cp_num$x.t
hist(time_data$arcsinh.rhizop_cp_num)

set.seed(13)
arcsinh.rhizop_cp_num <- arcsinh_x(growth_data$rhizop_cp_num)
growth_data$arcsinh.rhizop_cp_num <- arcsinh.rhizop_cp_num$x.t
hist(growth_data$arcsinh.rhizop_cp_num)



#__________ POSTHOC: Tukey HSD using grouped data _____________#

# Rhizosphere - 16S rRNA copies posthoc

# 1) Time Series 
rhizos_pwc_time <- time_data %>%
  group_by(age_days) %>%
  tukey_hsd(log.rhizos_cp_num ~ group)
rhizos_pwc_time
# 7        group Control Zeb             0  -0.724    -1.32    -0.132  0.0225 *           
# 14       group Control Zeb             0  -0.279    -0.540   -0.0189 0.0385 *

# 2) Growth Series 
rhizos_pwc_growth <- growth_data %>%
  group_by(stage) %>%
  tukey_hsd(arcsinh.rhizos_cp_num ~ group)
rhizos_pwc_growth
# V1    group Control Zeb             0    2.11     0.437    3.79   0.0197 *           



# Rhizoplane - 16S rRNA copies posthoc

# 1) Time Series 
rhizop_pwc_time <- time_data %>%
  group_by(age_days) %>%
  tukey_hsd(arcsinh.rhizop_cp_num ~ group)
rhizop_pwc_time
# 49       group Control Zeb             0   -1.24   -2.18      -0.301 0.016  *           
# 63       group Control Zeb             0    1.47    0.0361     2.90  0.0457 *  

# 2) Growth Series 
rhizop_pwc_growth <- growth_data %>%
  group_by(stage) %>%
  tukey_hsd(arcsinh.rhizop_cp_num ~ group)
rhizop_pwc_growth #NS

rhizop_pwc_growth_wil_test <- growth_data %>%
  group_by(stage) %>%
  wilcox_test(rhizop_cp_num ~ group)
rhizop_pwc_growth_wil_test #NS




# MAKE PLOTS FOR 16S rRNA GENE COPY NUMBER

# load the summarized data set
data_rhizo_assembly.tdy.df <- read.csv("/Users/emiliedehon/Nextcloud/Documents/PROJECT/Rhizosphere_assembly_Common_Bean/R_Analysis_Files/data_rhizo_assembly.tdy.df.csv", sep=",", row.names = 1)


# define the x-axis time scales for the control and delay treatment
ctr.t <- c(" ","3","7"," ","14", "21"," "," "," ", "35"," "," ", "49"," ", "63")
zeb.t <- c("3","7","14","21"," ", "35"," "," ","49"," "," ", "63"," ", " "," ")


# line chart - separated
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



# 1. Rhizosphere - 16S rRNA gene copies with standard error

# Rhizosphere 16S copies stats annotation for growth series
rhizos.cp_anno_LC_growth <- data.frame(xstar = c(3.93, 4.7, 5.7, 6.72, 8, 11.2, 14.2), 
                                       ystar = c(0, 0, 0, 0, 0, 0, 0), 
                                       lab = c( "*", "ns", "ns", "ns", "ns", "ns", "ns"))

# Rhizosphere 16S copies stats annotation for time series
rhizos.cp_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4, 6, 9, 12), x2 = c(2, 3, 5, 6, 10, 13, 15), 
                                     y1 = c(6, 7.7, 8.3, 9.1, 10, 11, 11.4), y2 = c(6.2, 7.9, 8.5, 9.3, 10.2, 11.2, 11.6), 
                                     xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(6.7, 8.2, 8.65, 9.9, 10.7, 11.8, 12.1),
                                     lab = c( "ns", "*", "*", "ns", "ns", "ns", "ns"))


# 16S rRNA copies line chart

data_rhizo_assembly.tdy.df$cp.growth.series <- "A1. Growth Series"

rhizos_cp_LC_growth <- ggplot(data_rhizo_assembly.tdy.df, 
                              aes(x=factor(stage_all_ed, level=stage_all_order), 
                                  y=Mean.cp.rhizos, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.cp.rhizos.low), ymax=ifelse(series == "time", NA, SE.cp.rhizos.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(0, 17))+
  ylab("16S rRNA Gene<br>(copy number g<sup>-1</sup> soil)")+
  labs(title="A. Rhizosphere", x = "Plant Growth Stage", linetype = "Series", color= "Treatment") +
  facet_grid(. ~ cp.growth.series)+
  scale_x_discrete(expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "solid", "time" = "blank")) +
  scale_shape_manual(values =c("growth" = 16, "time" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_markdown(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position = c(x=0.13, y=0.8),
        legend.box.background = element_rect(color = "black", size = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.8,xmax=-0.2,ymin=-9.8,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-9.8,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-9.8,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-1.5,xmax=-0.2,ymin=-12.2,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-12.2,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-12.2,ymax=-0.01)+
  annotation_custom(anno.title,xmin=-3,xmax=-2.5,ymin=-11,ymax=-0.01)+
  theme(plot.margin = unit(c(1,1,1,2), "lines")) +
  geom_text(inherit.aes=FALSE, data = rhizos.cp_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")

rhizos_cp_LC_growth

# Time Series
data_rhizo_assembly.tdy.df$cp.time.series <- "A2. Time Series"

rhizos_cp_LC_time <- ggplot(data_rhizo_assembly.tdy.df, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.cp.rhizos, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.cp.rhizos.low), ymax=ifelse(series == "growth", NA, SE.cp.rhizos.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6) +
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+
  theme_bw() +
  scale_y_continuous(limits = c(0, 17))+
  ylab("16S rRNA Gene<br>(copy number g<sup>-1</sup> soil)")+
  labs(x = "Plant Growth Stage", linetype = "Series", color= "Treatment") +
  facet_grid(. ~ cp.time.series)+
  scale_x_discrete(labels = c("Pre_emergence" = "Pre\nemergence"),expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "blank", "time" = "dashed")) +
  scale_shape_manual(values =c("time" = 17, "growth" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_markdown(size=20),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype=guide_legend(order=2), 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1,xmax=-0.2,ymin=-10.8,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-10.8,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-10.8,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-0.8,xmax=-0.2,ymin=-13.7,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-13.7,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-13.7,ymax=-0.01)+
  annotation_custom(anno.title,xmin=-2,xmax=-2.5,ymin=-12,ymax=-0.01)+ 
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2), "lines")) + 
  geom_text(inherit.aes=FALSE, data = rhizos.cp_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  geom_segment(inherit.aes=FALSE, data =rhizos.cp_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rhizos.cp_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rhizos.cp_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

rhizos_cp_LC_time

# combine figures
library(patchwork)

Copies.combined <- rhizos_cp_LC_growth / rhizos_cp_LC_time  #plot_layout(guides = "collect")
Copies.combined


# 2. Rhizoplane - 16S rRNA gene copies with standard error


# Rhizoplane 16S copies stat annotation for growth series
rhizop.cp_anno_LC_growth <- data.frame(xstar = c(3.9, 4.8, 5.8, 6.8, 8, 11.2, 14.2), 
                                       ystar = c(0, 0, 0, 0, 0, 0, 0), 
                                       lab = c( "ns", "ns", "ns", "ns", "ns", "ns", "ns"))

# Rhizoplane 16S copies stat annotation for time series
rhizop.cp_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4, 6, 9, 12), x2 = c(2, 3, 5, 6, 10, 13, 15), 
                                     y1 = c(5, 5.7, 8.5, 10.2, 11.5, 15, 16), y2 = c(5.2, 5.9, 8.7, 10.4, 11.7, 15.2, 16.2), 
                                     xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(5.8, 6.5, 9.3, 11, 12.5, 15.4, 16.4),
                                     lab = c( "ns", "ns", "ns", "ns", "ns", "*", "*"))


# Rhizoplane 16S rRNA copies line chart
data_rhizo_assembly.tdy.df$cp.rhizop.growth.series <- "B1. Growth Series"

rhizop_cp_LC_growth <- ggplot(data_rhizo_assembly.tdy.df, 
                              aes(x=factor(stage_all_ed, level=stage_all_order), 
                                  y=Mean.cp.rhizop, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.cp.rhizop.low), ymax=ifelse(series == "time", NA, SE.cp.rhizop.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(0, 17))+
  ylab("16S rRNA Gene (copy number g<sup>-1</sup> soil)")+
  labs(title="B. Rhizoplane", x = "Plant Growth Stage", linetype = "Series", color= "Treatment") +
  facet_grid(. ~ cp.rhizop.growth.series)+
  scale_x_discrete(expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "solid", "time" = "blank")) +
  scale_shape_manual(values =c("growth" = 16, "time" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_blank(),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        legend.box.background = element_rect(color = "black", size = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.8,xmax=-0.2,ymin=-9.8,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-9.8,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-9.8,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-9.8,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-1.5,xmax=-0.2,ymin=-12.2,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-12.2,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-12.2,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-12.2,ymax=-0.01)+
  annotation_custom(anno.title,xmin=-3,xmax=-2.5,ymin=-11,ymax=-0.01)+ 
  theme(plot.margin = unit(c(1,1,1,2), "lines"))+ #1,2,4.6,2.6
  geom_text(inherit.aes=FALSE, data = rhizop.cp_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")


rhizop_cp_LC_growth


# combine figures

# Time Series
data_rhizo_assembly.tdy.df$cp.rhizop.time.series <- "B2. Time Series"

rhizop_cp_LC_time <- ggplot(data_rhizo_assembly.tdy.df, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.cp.rhizop, colour = group, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.cp.rhizop.low), ymax=ifelse(series == "growth", NA, SE.cp.rhizop.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group.series,linetype =series),linewidth=0.6) +
  geom_point(aes(group=group.series,col = group), size=3, alpha=0.4)+
  theme_bw() +
  scale_y_continuous(limits = c(0, 17))+
  ylab("16S rRNA Gene (copy number g<sup>-1</sup> soil)")+
  labs(x = "Plant Growth Stage", linetype = "Series", color= "Treatment") +
  facet_grid(. ~ cp.rhizop.time.series)+
  scale_x_discrete(labels = c("Pre_emergence" = "Pre\nemergence"),expand = c(-1.136, 1.136))+
  scale_color_manual(values = group.series.col, labels = c("Control", "Delayed"))+
  scale_linetype_manual(values = c("growth" = "blank", "time" = "dashed")) +
  scale_shape_manual(values =c("time" = 17, "growth" = NA))+
  theme(plot.title = element_text(size=25, face = "bold"),
        strip.text.x = element_text(size = 24),
        axis.text=element_text(size=17),
        axis.title.y =element_blank(),
        axis.title.x = element_text(size=20, vjust = -0.75),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15),
        legend.position="none",
        panel.border = element_rect(size = 0.3))+
  guides(color = guide_legend(order=1),
         linetype=guide_legend(order=2), 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1,xmax=-0.2,ymin=-10.8,ymax=-0.01) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-10.8,ymax=-0.01) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-10.8,ymax=-0.01)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-10.8,ymax=-0.01)+
  annotation_custom(zeb.day,xmin=-0.8,xmax=-0.2,ymin=-13.7,ymax=-0.01) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-13.7,ymax=-0.01) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-13.7,ymax=-0.01)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-13.7,ymax=-0.01)+
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2), "lines")) + 
  geom_text(inherit.aes=FALSE, data = rhizop.cp_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="grey38") +
  geom_segment(inherit.aes=FALSE, data =rhizop.cp_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rhizop.cp_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rhizop.cp_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

rhizop_cp_LC_time

# combine figures
Copies.combined.rhizop <- rhizop_cp_LC_growth / rhizop_cp_LC_time  
Copies.combined.rhizop

# combines rhizosphere and rhizoplane together
cp.all <- Copies.combined | Copies.combined.rhizop

# save figure
setwd('/Users/emiliedehon/Nextcloud/Documents/PROJECT/Rhizosphere_assembly_Common_Bean/Figures/')
ggsave("CPall.tiff",
       cp.all , device = "tiff",
       width = 16.5, height =12, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")

