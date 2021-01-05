##################################################
#Environmental parameters: Dataset 2020
##################################################

#Load packages
library(ggplot2)
library(rstatix);packageVersion("rstatix")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


#############
#Import data
#############

#Import Sample Data - "metadata"
all.data <- read.csv("./data/Metadata_2015_2020.csv", h=T, sep = ",")

# Checking data structure
summary(all.data)
str(all.data)

#Edit metadeta dataframe
all.data <- filter(all.data, Dataset == "2020")
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("bottom", "surface"))

#Make separated files for surface 
all.data.surf=filter(all.data, Depth == "surface")
all.data.bot=filter(all.data, Depth == "bottom")


#############
#Plots
#############

#Graph A: Temperature and dissolved oxygen
all.data.a <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = Temperature_sea, fill = Location), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_point(mapping = aes(x = Season, y = Dissolved_Oxygen*4, shape = "Dissolved oxygen", group = Location), position = position_dodge(width=.5), stat = 'identity', size=1, show.legend = FALSE) +
  scale_y_continuous(name = "Temperature [°C]", limits = c(0,35), 
                     sec.axis = sec_axis(~./4, name = expression (paste("DO [mg ", L^-1, "]")))) +
  scale_fill_manual(values=location.col)  +
  theme_bw() +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,0,0,0), "lines"),legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9)) +
  facet_wrap(vars(Depth)) +
  labs(shape = NULL, fill = NULL)
all.data.a 

#Graph B: Salinity
all.data.b <- ggplot(data = all.data) +
  geom_point(mapping = aes(x = Season, y = Salinity, colour = Location), position = position_dodge(width=.5), stat = 'identity', show.legend = FALSE) +
  scale_y_continuous(name = "Salinity", limits = c(30,40)) +
  scale_colour_manual(values=location.col)   +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0,0,0), "lines"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9)) +
  labs(shape = NULL)
all.data.b 

#Graph C: DOC (y) / C/N (y2)
all.data.c <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = DOC, fill = Location), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  geom_point(mapping = aes(x = Season, y = C_N*9, shape = "C:N", group = Location), position = position_dodge(width=.5), stat = 'identity', size =1, show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("DOC [µmol ", L^-1, "]")), limits = c(0,200),
                     sec.axis = sec_axis(~./9, name = expression (paste("C:N")))) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,0,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.c 

#Graph D: TDN
all.data.d <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = TDN, fill = Location), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("TDN [µmol ", L^-1, "]")), limits = c(0,30)) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0,0,0), "lines"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.d 

#Graph E: NH4 (y1) / NO2+NO3 (y2)
all.data.e <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = NH4, fill = Location), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
  geom_point(aes(x = Season, y = (NO2+NO3)/2, shape = "NO2+NO3", group = Location), position = position_dodge(width=.5), stat = 'identity', size = 1, show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("NH " [4]^"  +", " [µmol ", L^-1, "]")), limits = c(0,4.7),
                     sec.axis = sec_axis(~.*2, name = expression (paste("NO " [2]^"  -", " + NO " [3]^"  -", " [µmol ", L^-1, "]")))) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,0,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.e

#Graph F: PO4 (y1) 
all.data.f <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = PO4, fill = Location), stat = 'identity', width = .5, position = position_dodge(), show.legend = FALSE) +
 # geom_point(aes(x = Season, y = SiO3/5, shape = "SiO3", group = Location), position = position_dodge(width=.5), stat = 'identity') +
  scale_y_continuous(name = expression (paste("PO"[4 ]^"  3+", " [µmol ", L^-1, "]")), limits = c(0,1)) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,0,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=9), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.f

#Graph G: Bacterial abundance
all.data.g <- ggplot(data = all.data, aes(fill=Location)) +
  geom_bar(mapping = aes(x = Season, y = BA_FC), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_errorbar(mapping = aes(x = Season, ymin = BA_FC - BA_FC_STDEV, ymax = BA_FC + BA_FC_STDEV), position = position_dodge(width = .5), stat = 'identity', width = 0.25) +
  scale_y_continuous(name = expression(paste("Bacterial abundance [N cell ", L^-1, "]"))) +
  scale_fill_manual(values=location.col)  +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0,0,0), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.g

#Graph H: Bacterial carbon production
all.data.h <- ggplot(data = all.data, aes(fill=Location)) +
  geom_bar(mapping = aes(x = Season, y = BCP_C), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_errorbar(mapping = aes(x = Season, ymin = BCP_C - BCP_C_STDEV, ymax = BCP_C + BCP_C_STDEV), position = position_dodge(width = .5), stat = 'identity', width = 0.25) +
  scale_y_continuous(name = expression(paste("BCP [µ C ", L^-1, " ", d^-1, "]"))) +
  scale_fill_manual(values=location.col)  +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "right", plot.margin = unit(c(0,0,0,0), "lines"), 
        axis.text.x = element_text(angle = 50, hjust = 1), axis.title = element_text(size=8), strip.text = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.h

ggarrange(all.data.a, all.data.b, all.data.c, all.data.d, all.data.e,  all.data.f, all.data.h, all.data.g,
          labels = c("B", "C", "D", "E", "F", "G", "H", "I"), font.label = list(size=9), nrow = 4, ncol = 2, align = c("v"))


##################################
#Supplementary figure: Chl a, SiO3
##################################

#Supplementary 1: SiO3
all.data.s1 <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = SiO3, fill = Location), stat = 'identity', width = .5, position = position_dodge()) +
  scale_y_continuous(name = expression (paste("SiO"[3 ], " [µmol ", L^-1, "]")), limits = c(0,15)) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  labs(shape = NULL, fill = NULL)
all.data.s1

#Supplementary 2:
all.data.s2 <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = Chl_A.1, fill = Location), stat = 'identity', width = .5, position = position_dodge()) +
  scale_y_continuous(name = expression (paste("Chl A [µg ", L^-1, "]")), limits = c(0,2)) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) +
  labs(shape = NULL, fill = NULL)
all.data.s2


##################################################
#Supplementary figure: PCA on environmental matrix
##################################################

#From https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples
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
envpar_corr.pvalues <- envpar_corr %>% cor_get_pval()

#plot
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)



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


# CLEAN UP #################################################

# Clear plots
dev.off() 

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
