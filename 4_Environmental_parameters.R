##################################################
#Environmental parameters: Dataset 2020
##################################################

#Load packages
library(ggplot2)
library(rstatix);packageVersion("rstatix")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
library(ggpubr)

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


#############
#Import data
#############

#Import Sample Data - "metadata"
all.data <- read.csv("./data/Metadata_2020.csv", sep = ",")


# Checking data structure
summary(all.data)
str(all.data)

#Edit metadeta dataframe
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("bottom", "surface"))

#Make separated files for surface 
all.data.surf=filter(all.data, Depth == "surface")
all.data.bot=filter(all.data, Depth == "bottom")


###################################
#Figure 1: Environmental parameters
###################################

#Graph Temperature
all.data.T <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Temperature_sea, fill = Season), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  scale_y_continuous(name = "Temperature [°C]", limits = c(0,30)) +
  scale_fill_manual(values=season.col)  +
  theme_bw() +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,2,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9),panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL) + 
  facet_wrap(vars(Depth))
all.data.T 

#Graph Dissolved oxygen
all.data.DO <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Dissolved_Oxygen, fill = Season), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("DO [mg ", L^-1, "]")), limits = c(0,10)) +
  scale_fill_manual(values=season.col)  +
  theme_bw() +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,2,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9),panel.grid = element_blank()) +
  facet_wrap(vars(Depth)) +
  labs(shape = NULL, fill = NULL)
all.data.DO

#Graph DOC
all.data.doc <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = DOC, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("DOC [µmol ", L^-1, "]")), limits = c(0,150)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,2,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9),panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.doc 

##Joined graphs for first line
graphs1 <- ggarrange(all.data.T, all.data.DO, all.data.doc,
                     labels = c("A", "B", "C"), font.label = list(size=9), nrow = 1, ncol = 3, align = c("hv"))
graphs1

#Graph NH4
all.data.nh4 <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = NH4, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("NH " [4]^"  +", " [µmol ", L^-1, "]")), limits = c(0,4)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank(), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.nh4

#Graph NO2+NO3
all.data.no <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = NO2+NO3, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("NO " [2]^"  -", " + NO " [3]^"  -", " [µmol ", L^-1, "]")), limits = c(0, 10)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank(), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.no

#Graph PO4
all.data.po4 <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = PO4, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("PO"[4 ]^"  3+", " [µmol ", L^-1, "]")), limits = c(0,1)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank(), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.po4

###Joined for inorganic nutrients
inorganic <- ggarrange(all.data.nh4, all.data.no, all.data.po4,
                       labels = c("D", "E", "F"), font.label = list(size=9), nrow = 1, ncol = 3, align = c("hv"))
inorganic

#Graph Bacterial abundance
all.data.ba <- ggplot(data = all.data, aes(fill=Season)) +
  geom_bar(mapping = aes(x = Location, y = BA_FC), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_errorbar(mapping = aes(x = Location, ymin = BA_FC - BA_FC_STDEV, ymax = BA_FC + BA_FC_STDEV), position = position_dodge(width = .5), stat = 'identity', width = 0.25) +
  scale_y_continuous(name = expression(paste("Bacterial abundance [N cell ", L^-1, "]"))) +
  scale_fill_manual(values=season.col)  +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,1,0.5,1), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.ba

#Graph Bacterial carbon production
all.data.bcp <- ggplot(data = all.data, aes(fill=Season)) +
  geom_bar(mapping = aes(x = Location, y = BCP_C), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_errorbar(mapping = aes(x = Location, ymin = BCP_C - BCP_C_STDEV, ymax = BCP_C + BCP_C_STDEV), position = position_dodge(width = .5), stat = 'identity', width = 0.25) +
  scale_y_continuous(name = expression(paste("BCP [µ C ", L^-1, " ", d^-1, "]"))) +
  scale_fill_manual(values=season.col)  +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,1,0.5,1), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.bcp

#Graph for bacteria
bacteria <- ggarrange(all.data.ba, all.data.bcp,
                      labels = c("G", "H"), font.label = list(size=9), nrow = 1, ncol = 2, align = c("hv"))
bacteria

figure1 <- ggarrange(graphs1, inorganic, bacteria, nrow = 3, ncol = 1, align = c("v"))

ggsave("./output_graphs/EnvironmentalParameters.pdf", figure1)


#######################
#Supplementary figure 1
#######################

#Supplementary 1:
#Graph Salinity
all.data.sal <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Salinity, fill = Season), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  scale_y_continuous(name = "Salinity", limits = c(0,40)) +
  scale_fill_manual(values=season.col)   +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0.5,1,2), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  facet_wrap(vars(Depth)) +
  labs(shape = NULL, fill = NULL)
all.data.sal 


#Graph TDN
all.data.tdn <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = TDN, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("TDN [µmol ", L^-1, "]")), limits = c(0,30)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0.5,1,2), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.tdn 

organic <- ggarrange(all.data.doc, all.data.tdn,
                     labels = c("A", "B"), font.label = list(size=9), nrow = 1, ncol = 2, align = c("hv"))
organic

#Graph C/N
all.data.cn <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = C_N, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = "C:N", limits = c(0,30)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0.5,1,2), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.cn 

#Graph SiO3
all.data.sio3 <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = SiO3, fill = Season), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("SiO"[3 ], " [µmol ", L^-1, "]")), limits = c(0,15)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0.5,1,2), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.sio3

#Graph Chla:
all.data.chla <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Chl_A.1, fill = Season), stat = 'identity', width = .5, position = position_dodge()) +
  scale_y_continuous(name = expression (paste("Chl A [µg ", L^-1, "]")), limits = c(0,2)) +
  scale_fill_manual(values=season.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "bottom", plot.margin = unit(c(0,0.5,1,2), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), panel.grid = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.chla

sup1 <- ggarrange(all.data.sal, all.data.tdn,all.data.cn, all.data.sio3, all.data.chla,
                  labels = c("A", "B", "C", "D", "E"), font.label = list(size=9), nrow = 3, ncol = 2, align = c("hv"))
sup1

ggsave("./output_graphs/EnvironmentalParameters_Sup.pdf", sup1)


###################################################
#Supplementary figure 2: Pearson correlation matrix
###################################################

#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import env.par. data
metadata.raw <- all.data
env.par <- c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "TDN", "C_N", "NO2_NO3", "PO4", "SiO3", "NH4", "BA_FC", "BCP_C", "Coliform_bacteria")

#scale all parameters
metadata.scaled.SRF <- metadata.raw %>% 
  mutate_at(env.par, scale_par)%>%
  as.data.frame()

#Select parameters
envpar_corr <- metadata.scaled.SRF %>% 
  select(env.par)%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr %>% cor_get_pval()

# Replacing correlation coefficients by symbols
envpar_corr  %>%   cor_as_symbols() %>%
  pull_lower_triangle()

#Mark significant correlations
signif_corr <- envpar_corr %>%
  cor_mark_significant()

#plot
sup2 <- envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)
sup2

ggsave("./output_graphs/CorrelationMatrix_Sup.pdf", sup2)


###########
#Statistic
###########

#Shapiro-Wilk's test
all.data %>% shapiro_test(Temperature_sea, Salinity, Dissolved_Oxygen, DOC, C_N, TDN, BCP_C, BA_FC, NO2_NO3, NO2, NO3, NH4, PO4)

#C_N, DOC - normal distribution - ANOVA
all.data %>% anova_test(DOC ~ Depth + Season + Location) 
all.data %>% tukey_hsd(DOC ~ Depth + Season + Location)

all.data %>% anova_test(C_N ~ Depth + Season + Location) 
all.data %>% tukey_hsd(C_N ~ Depth + Season + Location)

#Other - Kruskal-Wallis test
lapply(all.data[,c("Temperature_sea", "Salinity", "Dissolved_Oxygen", "TDN", "NO2_NO3", "NH4", "PO4")], function(x) kruskal.test(x ~ all.data$Season))
lapply(all.data[,c("Temperature_sea", "Salinity", "Dissolved_Oxygen", "TDN", "NO2_NO3", "NH4", "PO4")], function(x) kruskal.test(x ~ all.data$Location))
lapply(all.data[,c("Temperature_sea", "Salinity", "Dissolved_Oxygen", "TDN", "NO2_NO3", "NH4", "PO4")], function(x) kruskal.test(x ~ all.data$Depth))

lapply(all.data[,c("BCP_C", "BA_FC")], function(x) kruskal.test(x ~ all.data$Season))
lapply(all.data[,c("BCP_C", "BA_FC")], function(x) kruskal.test(x ~ all.data$Location))
lapply(all.data[,c("BCP_C", "BA_FC")], function(x) kruskal.test(x ~ all.data$Depth))

##Post-hoc:
lapply(all.data[,c("BCP_C", "BA_FC")], function(x) pairwise.wilcox.test(x, all.data$Season, p.adjust.method = "BH"))
lapply(all.data[,c("BCP_C", "BA_FC")], function(x) pairwise.wilcox.test(x, all.data$Location, p.adjust.method = "BH"))
lapply(all.data[,c("BCP_C", "BA_FC")], function(x) pairwise.wilcox.test(x, all.data$Depth, p.adjust.method = "BH"))


# CLEAN UP #################################################

# Clear plots
dev.off() 

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
