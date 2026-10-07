#title: "Rhizo Assembly Alpha Diversity Analysis"
#author: "Ludivine Guigard, based on Fina Bintarti"
#date: "05-10-2026"

  #Loading packages
library(multcomp)
library(car)
library(rstatix)
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
library(lme4)
library(nlme)
library(DHARMa)
library(phyloseq)

  #16S
getwd()
# load the multi-rarefied phyloseq object from Sam Barnett
RA_phyloseq_multrare <- readRDS("RA_phyloseq_multirarefied_Sam_16S.rds")
RA_phyloseq_multrare
# Make data frame of the ASV table
multrare.asv.df <-  as.data.frame(otu_table(RA_phyloseq_multrare))
str(multrare.asv.df)
# Richness
s <- specnumber(multrare.asv.df, MARGIN = 2) # richness
rich.df <- as.data.frame(s) 
# Read the updated metadata
metadata_update_edit <- read.csv("metadata_update_edit.csv", row.names = 1)
# Check if row names are identical (ignoring order)
setequal(rownames(metadata_update_edit), rownames(rich.df)) #TRUE
# Set SampleID
rich.df <- rownames_to_column(rich.df, var = "SampleID")
# Add column 's' (richness data) to the metadata_update_edit
metadata_update_edit$s <- rich.df$s[match(metadata_update_edit$SampleID, rich.df$SampleID)]
# Change necessary variables into factor
str(metadata_update_edit)
metadata_update_edit$s <- as.numeric(metadata_update_edit$s)
metadata_update_edit[sapply(metadata_update_edit, is.integer)] <- lapply(metadata_update_edit[sapply(metadata_update_edit, is.integer)], as.factor)
metadata_update_edit[sapply(metadata_update_edit, is.character)] <- lapply(metadata_update_edit[sapply(metadata_update_edit, is.character)], as.factor)
str(metadata_update_edit)
# Make a new data frame 
metadata_alphadiv <- metadata_update_edit
###############################################################################
## STATS
# Subset compartments
rsphere.df <- metadata_alphadiv %>% 
  filter(compartment == "rhizosphere")
rplane.df <- metadata_alphadiv %>% 
  filter(compartment == "rhizoplane")


# I. ANOVA - Alpha Diversity: Rhizosphere

# 1.) Time Series

# between group/treatment and age (days) in time series

# Keep only Time Series
rs_time <- rsphere.df %>% 
  filter(series == "time")
rs_time
rs_time <- droplevels(rs_time)
str(rs_time)

rs.rich.time.aov <- aov(s ~ treatment*age_days, data=rs_time)
summary(rs.rich.time.aov) #NS
#                    Df  Sum Sq Mean Sq F value Pr(>F)
#treatment           1   26880   26880   0.165  0.686
#age_days            6  754282  125714   0.773  0.595
#treatment:age_days  6  750128  125021   0.768  0.598

# Normality
hist(residuals(rs.rich.time.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rs.rich.time.aov))
qqline(residuals(rs.rich.time.aov), col = "red")
shapiro.test(residuals(rs.rich.time.aov)) # fine
# Homoscedasticity
plot(fitted(rs.rich.time.aov), residuals(rs.rich.time.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rs_time) #fine
leveneTest(s ~ age_days, data = rs_time) #fine
# POST-HOC ANALYSIS
# Tukey HSD
rs_rich_time.pwc <- rs_time %>%
  group_by(age_days) %>%
  tukey_hsd(s ~ treatment)
rs_rich_time.pwc #NS

# 2.) Growth Series

# between group/treatment and growth stage in growth series

# Keep only Growth Series
rs_growth <- rsphere.df %>% 
  filter(series == "growth")
rs_growth
rs_growth <- droplevels(rs_growth)
str(rs_growth)

rs.rich.growth.aov <- aov(s ~ treatment*stage, data=rs_growth)
summary(rs.rich.growth.aov) #NS
#.                Df  Sum Sq Mean Sq F value Pr(>F)
#treatment        1    8970    8970   0.063  0.803
#stage            6  720295  120049   0.838  0.546
#treatment:stage  6  527302   87884   0.613  0.719

# Normality
hist(residuals(rs.rich.growth.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rs.rich.growth.aov))
qqline(residuals(rs.rich.growth.aov), col = "red")
shapiro.test(residuals(rs.rich.growth.aov)) # fine
# Homoscedasticity
plot(fitted(rs.rich.growth.aov), residuals(rs.rich.growth.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rs_growth) #fine
leveneTest(s ~ stage, data = rs_growth) #fine

# POST-HOC ANALYSIS
# tukey hsd
rs_rich_growth.pwc <- rs_growth %>%
  group_by(stage) %>%
  tukey_hsd(s ~ treatment)
rs_rich_growth.pwc #NS


# II. ANOVA - Alpha Diversity: Rhizoplane

# 1.) Time Series

# between group/treatment and age (days) in time series

# Keep only Time Series
rp_time <- rplane.df %>% 
  filter(series == "time")
rp_time
rp_time <- droplevels(rp_time)
str(rp_time)

rp.rich.time.aov <- aov(s ~ treatment*age_days, data=rp_time)
summary(rp.rich.time.aov) 
#                    Df  Sum Sq Mean Sq F value  Pr(>F)   
#treatment           1  126133  126133   1.378 0.24567   
#age_days            6 1781994  296999   3.245 0.00863 **
#treatment:age_days  6 1773320  295553   3.229 0.00888 **

# Normality
hist(residuals(rp.rich.time.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rp.rich.time.aov))
qqline(residuals(rp.rich.time.aov), col = "red")
shapiro.test(residuals(rp.rich.time.aov)) # fine
# Homoscedasticity
plot(fitted(rp.rich.time.aov), residuals(rp.rich.time.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rp_time) #fine
leveneTest(s ~ age_days, data = rp_time) #fine
# POST-HOC ANALYSIS
# Tukey HSD
rp_rich_time.pwc <- rp_time %>%
  group_by(age_days) %>%
  tukey_hsd(s ~ treatment)
rp_rich_time.pwc 
# 35       treatment Control Zeb             0     725.     121.    1330.  0.0252 * 

# 2.) Growth Series

# between group/treatment and growth stage in growth series

# Keep only Growth Series
rp_growth <- rplane.df %>% 
  filter(series == "growth")
rp_growth
rp_growth <- droplevels(rp_growth)
str(rp_growth)

rp.rich.growth.aov <- aov(s ~ treatment*stage, data=rp_growth)
summary(rp.rich.growth.aov) 
#                 Df  Sum Sq Mean Sq F value  Pr(>F)   
#treatment        1   96582   96582   0.879 0.35294   
#stage            6 2863165  477194   4.341 0.00127 **
#treatment:stage  6  298850   49808   0.453 0.83952 

# Normality
hist(residuals(rp.rich.growth.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rp.rich.growth.aov))
qqline(residuals(rp.rich.growth.aov), col = "red")
shapiro.test(residuals(rp.rich.growth.aov)) # fine
# Homoscedasticity
plot(fitted(rp.rich.growth.aov), residuals(rp.rich.growth.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rp_growth) #fine
leveneTest(s ~ stage, data = rp_growth) #fine

# POST-HOC ANALYSIS
# tukey hsd
rp_rich_growth.pwc <- rp_growth %>%
  group_by(stage) %>%
  tukey_hsd(s ~ treatment)
rp_rich_growth.pwc #NS

## ALPHA DIVERSITY PLOTS

# 1. Rhizosphere

# Load the rhizosphere dataframe
view(rsphere.df)
rsphere.df <- droplevels(rsphere.df)
str(rsphere.df)

# Tidy up the data and summarize
rsphere.df.tdy <- rsphere.df %>%
  group_by(treatment, series, group,stage,stage_all_ed,age_days) %>%
  dplyr::summarize(N = n(),
            Mean.rich=if(all(is.na(s))) NA_real_ else mean(s, na.rm = T),
            SD.rich=if(all(is.na(s))) NA_real_ else sd(s, na.rm = T),
            SE.rich.low = Mean.rich - (SD.rich/sqrt(N)),
            SE.rich.high = Mean.rich + (SD.rich/sqrt(N)))
str(rsphere.df.tdy)
# Re-order the growth stage
rsphere.df.tdy$stage_all_ed <- factor(rsphere.df.tdy$stage_all_ed, levels = c("PE","VE","VC",
                                                                                    "V1","V2","V3","V4","R1","R2",
                                                                                    "R3","R4","R5","R6","R7","RH"))
rsphere.df.tdy$stage_all_ed
stage_all_order <- c("PE", "VE","VC","V1","V2","V3", "V4", 
                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "RH")

# define the x-axis time scales for the control and delay treatment
ctr.t <- c(" ","3","7"," ","14", "21"," "," "," ", "35"," "," ", "49"," ", "63")
zeb.t <- c("3","7","14","21"," ", "35"," "," ","49"," "," ", "63"," ", " "," ")

# assign the x-axis color to match the legend
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
                   gp = gpar(lty=2,lwd = 1, col = "black", fill="#00000000"))
# border title
anno.title <- textGrob("Time\n(days)", gp=gpar(fontsize=14), rot=90)
# assign colors
group.series.col <- c("#0072B2","#D55E00")


# 1.) Rhizosphere - Growth Series

# Richness stat annotation for growth series
rs.rich_anno_LC_growth <- data.frame(xstar = c(3.72, 4.85, 5.85, 6.72, 8, 11.2, 14.2), 
                                     ystar = c(500, 500, 500, 500, 500, 500, 500), 
                                     lab = c( "ns", "ns", "ns", "ns", "ns", "ns", "ns"))
rs.rich_anno_LC_growth

rsphere.df.tdy$rich.growth.series <- "A1. Growth Series"

rs.RICH_LC_growth <- ggplot(rsphere.df.tdy, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.rich.low), ymax=ifelse(series == "time", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ rich.growth.series)+
  scale_y_continuous(limits = c(500, 2000))+
  ylab("Number of Observed ASVs")+
  labs(title="A. Rhizosphere", x = "Plant Growth Stage", color= "Treatment") +
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))
  #geom_text(inherit.aes=FALSE, data = rs.rich_anno_LC_growth, 
            #aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")
rs.RICH_LC_growth

# 2.) Rhizosphere - Time Series

# Richness stat annotation for time series
rs.rich_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4,6,9,12), x2 = c(2, 3, 5, 6,10,13,15), 
                                   y1 = c(1680, 1750, 1780, 1880, 1780, 1880, 1920), y2 = c(1700, 1770, 1800, 1900, 1800, 1900, 1940), 
                                   xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(1750, 1820, 1850, 1950, 1850 , 1980, 2000),
                                   lab = c( "ns", "ns", "ns", "ns", "ns", "ns", "ns"))
rs.rich_anno_LC_time

rsphere.df.tdy$rich.time.series <- "A2. Time Series"

rs.RICH_LC_time <- ggplot(rsphere.df.tdy, 
                          aes(x=factor(stage_all_ed, level=stage_all_order), 
                              y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.rich.low), ymax=ifelse(series == "growth", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(500, 2000))+
  facet_grid(. ~ rich.time.series)+
  ylab("Number of Observed ASVs")+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.2,ymin=-40,ymax=120) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=-40,ymax=120) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=-40,ymax=120)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=-40,ymax=120)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=-40,ymax=120)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=-40,ymax=120)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=-40,ymax=120)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=-40,ymax=120)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=-100,ymax=-100) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-100,ymax=-100) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-100,ymax=-100)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-100,ymax=-100)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-100,ymax=-100)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-100,ymax=-100)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-100,ymax=-100)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-100,ymax=-100)+
  annotation_custom(anno.title,xmin=-2.4,xmax=-2,ymin=-10,ymax=-100)+
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) 
  #geom_text(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  #geom_segment(inherit.aes=FALSE, data =rs.rich_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               #colour = "grey") +
  #geom_segment(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               #colour = "grey") +
  #geom_segment(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               #colour = "grey")

rs.RICH_LC_time

library("patchwork")
rs.combined <- rs.RICH_LC_growth / rs.RICH_LC_time  


# 2. Rhizoplane

# Load the rhizosphere dataframe
view(rplane.df)
rplane.df <- droplevels(rplane.df)
str(rplane.df)

# Tidy up the data and summarize
rplane.df.tdy <- rplane.df %>%
  group_by(treatment, series, group,stage,stage_all_ed,age_days) %>%
  dplyr::summarize(N = n(),
            Mean.rich=if(all(is.na(s))) NA_real_ else mean(s, na.rm = T),
            SD.rich=if(all(is.na(s))) NA_real_ else sd(s, na.rm = T),
            SE.rich.low = Mean.rich - (SD.rich/sqrt(N)),
            SE.rich.high = Mean.rich + (SD.rich/sqrt(N)))
str(rplane.df.tdy)
# Re-order the growth stage
rplane.df.tdy$stage_all_ed <- factor(rplane.df.tdy$stage_all_ed, levels = c("PE","VE","VC",
                                                                              "V1","V2","V3","V4","R1","R2",
                                                                              "R3","R4","R5","R6","R7","RH"))
rplane.df.tdy$stage_all_ed
stage_all_order <- c("PE", "VE","VC","V1","V2","V3", "V4", 
                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "RH")

# 1.) Rhizoplane - Growth Series

# Dry weight stat annotation for growth series
rp.rich_anno_LC_growth <- data.frame(xstar = c(3.72, 4.85, 5.85, 6.72, 8, 11.2, 14.2), 
                                     ystar = c(500, 500, 500, 500, 500, 500, 500), 
                                     lab = c( "", "", "", "", "", "", ""))
                                     #lab = c( "ns", "ns", "ns", "ns", "ns", "*", "ns"))

rplane.df.tdy$rich.growth.series <- "B1. Growth Series"

rp.RICH_LC_growth <- ggplot(rplane.df.tdy, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.rich.low), ymax=ifelse(series == "time", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ rich.growth.series)+
  scale_y_continuous(limits = c(500, 2000))+
  labs(title="B. Rhizoplane", x = "Plant Growth Stage", color= "Treatment") +
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
        legend.position = "none",
        legend.box.background = element_rect(color = "black", linewidth = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))+
  geom_text(inherit.aes=FALSE, data = rp.rich_anno_LC_growth, 
            aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")

rp.RICH_LC_growth


# 2.) Rhizoplane - Time Series

# Dry weight stat annotation for time series
rp.rich_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4,6,9,12), x2 = c(2, 3, 5, 6,10,13,15), 
                                   y1 = c(1680, 1750, 1790, 1880, 1780, 1880, 1920), y2 = c(1700, 1770, 1810, 1900, 1800, 1900, 1940), 
                                   xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(1750, 1820, 1830, 1950, 1850 , 1980, 2000),
                                   lab = c( "", "", "", "", "*", "", ""))
                                   #lab = c( "ns", "ns", "*", "ns", "*", "ns", "ns"))
rp.rich_anno_LC_time

rp.rich_anno_LC_time <- data.frame(x1 = c(6), x2 = c(10), 
                                   y1 = c(1840), y2 = c(1860), 
                                   xstar = c(8), ystar = c(1865),
                                   lab = c("*"))



rplane.df.tdy$rich.time.series <- "B2. Time Series"

rp.RICH_LC_time <- ggplot(rplane.df.tdy, 
                          aes(x=factor(stage_all_ed, level=stage_all_order), 
                              y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.rich.low), ymax=ifelse(series == "growth", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous(limits = c(500, 2000))+
  facet_grid(. ~ rich.time.series)+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.2,ymin=30,ymax=30) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=30,ymax=30) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=30,ymax=30)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=30,ymax=30)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=30,ymax=30)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=30,ymax=30)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=30,ymax=30)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=30,ymax=30)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=-120,ymax=-90) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=-120,ymax=-90) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=-120,ymax=-90)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=-120,ymax=-90)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=-120,ymax=-90)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=-120,ymax=-90)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=-120,ymax=-90)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=-120,ymax=-90)+
  #annotation_custom(anno.title,xmin=-2.4,xmax=-2,ymin=-250,ymax=-80)+
  annotation_custom(border) + coord_cartesian(clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) +
  geom_text(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
  geom_segment(inherit.aes=FALSE, data =rp.rich_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
               colour = "grey") +
  geom_segment(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "grey")

rp.RICH_LC_time

# Combine plots
rp.combined <- rp.RICH_LC_growth / rp.RICH_LC_time 

# Combine Rhizosphere and Rhizoplane Plots
RICH.all.16S <- rs.combined | rp.combined
RICH.all.16S
ggsave("16S.tiff",
       RICH.all.16S, device = "tiff",
       width = 16.5, height =12, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")


  #ITS2
  #Same base until the end where we combine figures so be careful: 
  #make sure you have the right figure parameters first

# load the multi-rarefied phyloseq object
getwd()
RA_phyloseq_multrare <- readRDS("RA_combined_final.rds")
RA_phyloseq_multrare
taxa_are_rows(RA_phyloseq_multrare)
#If taxa are rows = FALSE: need to change the richness calculation to MARGIN = 1
# Make data frame of the ASV table
multrare.asv.df <-  as.data.frame(otu_table(RA_phyloseq_multrare))
str(multrare.asv.df)
# Richness
s <- specnumber(multrare.asv.df, MARGIN = 1) # richness
rich.df <- as.data.frame(s) 
# Read the updated metadata
metadata_update_edit <- read.csv("metadata_ITS.csv", sep = ";", row.names = 1)
# Check if row names are identical (ignoring order)
setequal(rownames(metadata_update_edit), rownames(rich.df)) #TRUE
# Set SampleID
rich.df <- rownames_to_column(rich.df, var = "SampleID")
# Add column 's' (richness data) to the metadata_update_edit
metadata_update_edit$s <- rich.df$s[match(metadata_update_edit$SampleID, rich.df$SampleID)]
# Change necessary variables into factor
str(metadata_update_edit)
metadata_update_edit$s <- as.numeric(metadata_update_edit$s)
metadata_update_edit[sapply(metadata_update_edit, is.integer)] <- lapply(metadata_update_edit[sapply(metadata_update_edit, is.integer)], as.factor)
metadata_update_edit[sapply(metadata_update_edit, is.character)] <- lapply(metadata_update_edit[sapply(metadata_update_edit, is.character)], as.factor)
str(metadata_update_edit)
# Make a new data frame 
metadata_alphadiv <- metadata_update_edit
###############################################################################
## STATS
# Subset compartments
rsphere.df <- metadata_alphadiv %>% 
  filter(compartment == "rhizosphere")
rplane.df <- metadata_alphadiv %>% 
  filter(compartment == "rhizoplane")

# I. ANOVA - Alpha Diversity: Rhizosphere

# 1.) Time Series

# between group/treatment and age (days) in time series

# Keep only Time Series
rs_time <- rsphere.df %>% 
  filter(series == "time")
rs_time
rs_time <- droplevels(rs_time)
str(rs_time)

rs.rich.time.aov <- aov(s ~ treatment*age_days, data=rs_time)
summary(rs.rich.time.aov) #NS
#Df Sum Sq Mean Sq F value Pr(>F)
#treatment           1      1       1   0.000  0.993
#age_days            6  35118    5853   0.500  0.806
#treatment:age_days  6 103817   17303   1.478  0.203
#Residuals          55 643836   11706 

# Normality
hist(residuals(rs.rich.time.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rs.rich.time.aov))
qqline(residuals(rs.rich.time.aov), col = "red")
shapiro.test(residuals(rs.rich.time.aov)) # fine
# Homoscedasticity
plot(fitted(rs.rich.time.aov), residuals(rs.rich.time.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rs_time) #fine
leveneTest(s ~ age_days, data = rs_time) #fine
# POST-HOC ANALYSIS
# Tukey HSD
rs_rich_time.pwc <- rs_time %>%
  group_by(age_days) %>%
  tukey_hsd(s ~ treatment)
rs_rich_time.pwc #NS

# 2.) Growth Series

# between group/treatment and growth stage in growth series

# Keep only Growth Series
rs_growth <- rsphere.df %>% 
  filter(series == "growth")
rs_growth
rs_growth <- droplevels(rs_growth)
str(rs_growth)

rs.rich.growth.aov <- aov(s ~ treatment*stage, data=rs_growth)
summary(rs.rich.growth.aov) #tendancies?

#without log transfo
#Df Sum Sq Mean Sq F value Pr(>F)  
#treatment        1  42159   42159   3.245 0.0777 .
#stage            6 171644   28607   2.202 0.0581 .
#treatment:stage  6  69783   11631   0.895 0.5058  
#Residuals       50 649626   12993 

#Normality
hist(residuals(rs.rich.growth.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rs.rich.growth.aov))
qqline(residuals(rs.rich.growth.aov), col = "red")
shapiro.test(residuals(rs.rich.growth.aov)) # fine
# Homoscedasticity
plot(fitted(rs.rich.growth.aov), residuals(rs.rich.growth.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rs_growth) #0.02, should transform

rs_growth$log_s <- log(rs_growth$s)
rs.rich.growth.aov <- aov(log_s ~ treatment*stage, data=rs_growth)
summary(rs.rich.growth.aov)
#with the log:
#                Df Sum Sq Mean Sq F value Pr(>F)  
#treatment        1 0.0590 0.05897   2.442 0.1244  
#stage            6 0.2868 0.04780   1.980 0.0862 .
#treatment:stage  6 0.1134 0.01890   0.783 0.5875  
#Residuals       50 1.2073 0.02415                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Log_Normality
hist(residuals(rs.rich.growth.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rs.rich.growth.aov))
qqline(residuals(rs.rich.growth.aov), col = "red")
shapiro.test(residuals(rs.rich.growth.aov)) # fine
#Log_Homoscedasticity
plot(fitted(rs.rich.growth.aov), residuals(rs.rich.growth.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(log_s ~ treatment, data = rs_growth)
leveneTest(log_s ~ stage, data = rs_growth) #fine

# POST-HOC ANALYSIS
# tukey hsd
rs_rich_growth.pwc <- rs_growth %>%
  group_by(stage) %>%
  tukey_hsd(log_s ~ treatment)
rs_rich_growth.pwc #NS

# II. ANOVA - Alpha Diversity: Rhizoplane

# 1.) Time Series

# between group/treatment and age (days) in time series

# Keep only Time Series
rp_time <- rplane.df %>% 
  filter(series == "time")
rp_time
rp_time <- droplevels(rp_time)
str(rp_time)

rp.rich.time.aov <- aov(s ~ treatment*age_days, data=rp_time)
summary(rp.rich.time.aov) 
#Df Sum Sq Mean Sq F value Pr(>F)  
#treatment           1  63400   63400   3.197 0.0820 .
#age_days            6 363101   60517   3.052 0.0158 *
#treatment:age_days  6 176239   29373   1.481 0.2115  
#Residuals          37 733679   19829 

# Normality
hist(residuals(rp.rich.time.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rp.rich.time.aov))
qqline(residuals(rp.rich.time.aov), col = "red")
shapiro.test(residuals(rp.rich.time.aov)) # fine
# Homoscedasticity
plot(fitted(rp.rich.time.aov), residuals(rp.rich.time.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rp_time) #fine
leveneTest(s ~ age_days, data = rp_time) #fine
# POST-HOC ANALYSIS
# Tukey HSD
rp_rich_time.pwc <- rp_time %>%
  group_by(age_days) %>%
  tukey_hsd(s ~ treatment)
rp_rich_time.pwc #NS

# 2.) Growth Series

# between group/treatment and growth stage in growth series

# Keep only Growth Series
rp_growth <- rplane.df %>% 
  filter(series == "growth")
rp_growth
rp_growth <- droplevels(rp_growth)
str(rp_growth)

rp.rich.growth.aov <- aov(s ~ treatment*stage, data=rp_growth)
summary(rp.rich.growth.aov) 
#Df  Sum Sq Mean Sq F value Pr(>F)  
#treatment        1   53807   53807   1.507 0.2260  
#stage            6  532491   88748   2.485 0.0367 *
#treatment:stage  6   52257    8709   0.244 0.9593 
#Residuals       45 1607130   35714

# Normality
hist(residuals(rp.rich.growth.aov), main = "Histogram of Residuals", xlab = "Residuals") # fine
qqnorm(residuals(rp.rich.growth.aov))
qqline(residuals(rp.rich.growth.aov), col = "red")
shapiro.test(residuals(rp.rich.growth.aov)) # fine
# Homoscedasticity
plot(fitted(rp.rich.growth.aov), residuals(rp.rich.growth.aov), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
leveneTest(s ~ treatment, data = rp_growth) #fine
leveneTest(s ~ stage, data = rp_growth) #fine

# POST-HOC ANALYSIS
# tukey hsd
rp_rich_growth.pwc <- rp_growth %>%
  group_by(stage) %>%
  tukey_hsd(s ~ treatment)
rp_rich_growth.pwc #NS

## ALPHA DIVERSITY PLOTS

# 1. Rhizosphere

# Load the rhizosphere dataframe
view(rsphere.df)
rsphere.df <- droplevels(rsphere.df)
str(rsphere.df)

# Tidy up the data and summarize
rsphere.df.tdy <- rsphere.df %>%
  group_by(treatment, series, group,stage,stage_all_ed,age_days) %>%
  dplyr::summarize(N = n(),
                   Mean.rich=if(all(is.na(s))) NA_real_ else mean(s, na.rm = T),
                   SD.rich=if(all(is.na(s))) NA_real_ else sd(s, na.rm = T),
                   SE.rich.low = Mean.rich - (SD.rich/sqrt(N)),
                   SE.rich.high = Mean.rich + (SD.rich/sqrt(N)))
str(rsphere.df.tdy)
# Re-order the growth stage
rsphere.df.tdy$stage_all_ed <- factor(rsphere.df.tdy$stage_all_ed, levels = c("PE","VE","VC",
                                                                              "V1","V2","V3","V4","R1","R2",
                                                                              "R3","R4","R5","R6","R7","RH"))
rsphere.df.tdy$stage_all_ed
stage_all_order <- c("PE", "VE","VC","V1","V2","V3", "V4", 
                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "RH")

# define the x-axis time scales for the control and delay treatment
ctr.t <- c(" ","3","7"," ","14", "21"," "," "," ", "35"," "," ", "49"," ", "63")
zeb.t <- c("3","7","14","21"," ", "35"," "," ","49"," "," ", "63"," ", " "," ")

# assign the x-axis color to match the legend
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
                   gp = gpar(lty=2,lwd = 1, col = "black", fill="#00000000"))
# border title
anno.title <- textGrob("Time\n(days)", gp=gpar(fontsize=14), rot=90)
# assign colors
group.series.col <- c("#0072B2","#D55E00")


# 1.) Rhizosphere - Growth Series

# Richness stat annotation for growth series
rs.rich_anno_LC_growth <- data.frame(xstar = c(3.72, 4.85, 5.85, 6.72, 8, 11.2, 14.2), 
                                     ystar = c(500, 500, 500, 500, 500, 500, 500), 
                                     lab = c( "ns", "ns", "ns", "ns", "ns", "ns", "ns"))
rs.rich_anno_LC_growth

rsphere.df.tdy$rich.growth.series <- "A1. Growth Series"

rs.RICH_LC_growth <- ggplot(rsphere.df.tdy, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.rich.low), ymax=ifelse(series == "time", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ rich.growth.series)+
  scale_y_continuous(limits = c(500, 1300))+
  ylab("Number of Observed ASVs")+
  labs(title="A. Rhizosphere", x = "Plant Growth Stage", color= "Treatment") +
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))
#geom_text(inherit.aes=FALSE, data = rs.rich_anno_LC_growth, 
#aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")
rs.RICH_LC_growth

# 2.) Rhizosphere - Time Series

# Richness stat annotation for time series
rs.rich_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4,6,9,12), x2 = c(2, 3, 5, 6,10,13,15), 
                                   y1 = c(1680, 1750, 1780, 1880, 1780, 1880, 1920), y2 = c(1700, 1770, 1800, 1900, 1800, 1900, 1940), 
                                   xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(1750, 1820, 1850, 1950, 1850 , 1980, 2000),
                                   lab = c( "ns", "ns", "ns", "ns", "ns", "ns", "ns"))
rs.rich_anno_LC_time

rsphere.df.tdy$rich.time.series <- "A2. Time Series"

rs.RICH_LC_time <- ggplot(rsphere.df.tdy, 
                          aes(x=factor(stage_all_ed, level=stage_all_order), 
                              y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.rich.low), ymax=ifelse(series == "growth", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous()+
  facet_grid(. ~ rich.time.series)+
  ylab("Number of Observed ASVs")+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.2,ymin=255,ymax=255) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=255,ymax=255) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=255,ymax=255)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=255,ymax=255)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=255,ymax=255)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=255,ymax=255)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=255,ymax=255)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=255,ymax=255)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=180,ymax=180) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=180,ymax=180) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=180,ymax=180)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=180,ymax=180)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=180,ymax=180)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=180,ymax=180)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=180,ymax=180)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=180,ymax=180)+
  annotation_custom(anno.title,xmin=-2.2,xmax=-2.2,ymin=205,ymax=205)+
  annotation_custom(border) + 
  coord_cartesian(ylim = c(500, 1300),
                  clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) 
#geom_text(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
#geom_segment(inherit.aes=FALSE, data =rs.rich_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
#colour = "grey") +
#geom_segment(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
#colour = "grey") +
#geom_segment(inherit.aes=FALSE, data = rs.rich_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
#colour = "grey")

rs.RICH_LC_time

rs.combined <- rs.RICH_LC_growth / rs.RICH_LC_time  


# 2. Rhizoplane

# Load the rhizosphere dataframe
view(rplane.df)
rplane.df <- droplevels(rplane.df)
str(rplane.df)

# Tidy up the data and summarize
rplane.df.tdy <- rplane.df %>%
  group_by(treatment, series, group,stage,stage_all_ed,age_days) %>%
  dplyr::summarize(N = n(),
                   Mean.rich=if(all(is.na(s))) NA_real_ else mean(s, na.rm = T),
                   SD.rich=if(all(is.na(s))) NA_real_ else sd(s, na.rm = T),
                   SE.rich.low = Mean.rich - (SD.rich/sqrt(N)),
                   SE.rich.high = Mean.rich + (SD.rich/sqrt(N)))
str(rplane.df.tdy)
# Re-order the growth stage
rplane.df.tdy$stage_all_ed <- factor(rplane.df.tdy$stage_all_ed, levels = c("PE","VE","VC",
                                                                            "V1","V2","V3","V4","R1","R2",
                                                                            "R3","R4","R5","R6","R7","RH"))
rplane.df.tdy$stage_all_ed
stage_all_order <- c("PE", "VE","VC","V1","V2","V3", "V4", 
                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "RH")

# 1.) Rhizoplane - Growth Series

# Dry weight stat annotation for growth series
rp.rich_anno_LC_growth <- data.frame(xstar = c(3.72, 4.85, 5.85, 6.72, 8, 11.2, 14.2), 
                                     ystar = c(500, 500, 500, 500, 500, 500, 500), 
                                     lab = c( "", "", "", "", "", "", ""))
#lab = c( "ns", "ns", "ns", "ns", "ns", "*", "ns"))

rplane.df.tdy$rich.growth.series <- "B1. Growth Series"

rp.RICH_LC_growth <- ggplot(rplane.df.tdy, 
                            aes(x=factor(stage_all_ed, level=stage_all_order), 
                                y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "time", NA, SE.rich.low), ymax=ifelse(series == "time", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3, inherit.aes = T)+ 
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  facet_grid(. ~ rich.growth.series)+
  scale_y_continuous(limits = c(500, 1300))+
  labs(title="B. Rhizoplane", x = "Plant Growth Stage", color= "Treatment") +
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
        legend.position = "none",
        legend.box.background = element_rect(color = "black", linewidth = .5),
        legend.key.size = unit(1, "cm"),
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  theme(plot.margin = unit(c(1,1,1,2.3), "lines"))
#+geom_text(inherit.aes=FALSE, data = rp.rich_anno_LC_growth, 
#aes(x = xstar,  y = ystar, label = lab), size=7,vjust="inward",hjust="inward")

rp.RICH_LC_growth


# 2.) Rhizoplane - Time Series

# Dry weight stat annotation for time series
#rp.rich_anno_LC_time <- data.frame(x1 = c(1, 2, 3, 4,6,9,12), x2 = c(2, 3, 5, 6,10,13,15), 
# y1 = c(1680, 1750, 1790, 1880, 1780, 1880, 1920), y2 = c(1700, 1770, 1810, 1900, 1800, 1900, 1940), 
#xstar = c(1.5, 2.5, 4, 5, 8, 11, 13.5), ystar = c(1750, 1820, 1830, 1950, 1850 , 1980, 2000),
#lab = c( "", "", "", "", "", "", ""))
#lab = c( "ns", "ns", "*", "ns", "*", "ns", "ns"))
#rp.rich_anno_LC_time


rplane.df.tdy$rich.time.series <- "B2. Time Series"

rp.RICH_LC_time <- ggplot(rplane.df.tdy, 
                          aes(x=factor(stage_all_ed, level=stage_all_order), 
                              y=Mean.rich, colour = treatment, shape=series)) +
  geom_errorbar(aes(ymin=ifelse(series == "growth", NA, SE.rich.low), ymax=ifelse(series == "growth", NA, SE.rich.high)), width=.2,
                position=position_dodge(0.05), alpha=0.3)+
  geom_line(aes(group=group,linetype =series),linewidth=0.6)+ 
  geom_point(aes(group=group,col = treatment), size=3, alpha=0.4)+ 
  theme_bw() +
  scale_y_continuous()+
  facet_grid(. ~ rich.time.series)+
  labs(x = "Plant Growth Stage", color= "Treatment") +
  scale_x_discrete(expand = c(-1.136, 1.136))+
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
        panel.border = element_rect(linewidth = 0.3))+
  guides(color = guide_legend(order=1),
         linetype="none", 
         shape="none")+
  annotation_custom(ctr.day,xmin=-1.2,xmax=-0.2,ymin=255,ymax=255) +
  annotation_custom(ctr.3,xmin=2,xmax=2,ymin=255,ymax=255) +
  annotation_custom(ctr.7,xmin=3,xmax=3,ymin=255,ymax=255)+
  annotation_custom(ctr.14,xmin=5,xmax=5,ymin=255,ymax=255)+
  annotation_custom(ctr.21,xmin=6,xmax=6,ymin=255,ymax=255)+
  annotation_custom(ctr.35,xmin=10,xmax=10,ymin=255,ymax=255)+
  annotation_custom(ctr.49,xmin=13,xmax=13,ymin=255,ymax=255)+
  annotation_custom(ctr.63,xmin=15,xmax=15,ymin=255,ymax=255)+
  annotation_custom(zeb.day,xmin=-0.9,xmax=-0.2,ymin=180,ymax=180) +
  annotation_custom(zeb.3,xmin=1,xmax=1,ymin=180,ymax=180) +
  annotation_custom(zeb.7,xmin=2,xmax=2,ymin=180,ymax=180)+
  annotation_custom(zeb.14,xmin=3,xmax=3,ymin=180,ymax=180)+
  annotation_custom(zeb.21,xmin=4,xmax=4,ymin=180,ymax=180)+
  annotation_custom(zeb.35,xmin=6,xmax=6,ymin=180,ymax=180)+
  annotation_custom(zeb.49,xmin=9,xmax=9,ymin=180,ymax=180)+
  annotation_custom(zeb.63,xmin=12,xmax=12,ymin=180,ymax=180)+
  #annotation_custom(anno.title,xmin=-2.4,xmax=-2,ymin=-250,ymax=-80)+
  annotation_custom(border) + coord_cartesian(ylim = c(500, 1300),
                                              clip = "off")+
  theme(plot.margin = unit(c(1,1,4.2,2.3), "lines")) #+
#geom_text(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = xstar,  y = ystar, label = lab), size=7, color="black") +
#geom_segment(inherit.aes=FALSE, data =rp.rich_anno_LC_time, aes(x = x1, xend = x1, y = y1, yend = y2),
#             colour = "grey") +
#geom_segment(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = x2, xend = x2, y = y1, yend = y2),
#             colour = "grey") +
#geom_segment(inherit.aes=FALSE, data = rp.rich_anno_LC_time, aes(x = x1, xend = x2, y = y2, yend = y2),
#             colour = "grey")

rp.RICH_LC_time

# Combine plots
rp.combined <- rp.RICH_LC_growth / rp.RICH_LC_time 

# Combine Rhizosphere and Rhizoplane Plots
RICH.all.ITS <- rs.combined | rp.combined
RICH.all.ITS

  #Super Figure

title_row <- function(txt) {
  wrap_elements(full = grid::textGrob(txt, x = 0, hjust = 0,
                                      gp = grid::gpar(fontface = "bold", fontsize = 24)))
}

RICH.all <- (title_row("1. Number of observed bacterial ASVs (16S)") /
               wrap_elements(full = RICH.all.16S) /
               title_row("2. Number of observed fungal ASVs (ITS2)") /
               wrap_elements(full = RICH.all.ITS)) +
  plot_layout(heights = c(0.05, 0.8, 0.07, 0.8))

ggsave("Alpha_Div.tiff",
       RICH.all, device = "tiff",
       width = 19, height = 27, 
       units= "in", dpi = 300,
       compression="lzw", bg= "white")
