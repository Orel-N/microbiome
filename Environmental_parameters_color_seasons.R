#MICROBIOME - Dynamics of water biotic and abiotic paremeters
#Neža Orel, neza.orel@nib.si
#Script for graphs of accompining physical, chemical and biological properties.
#Data was collected during in situ survey 2020, project MIKROBIOM

############################################################################

#Set working directory
setwd("C:/Users/nezao/Documents/5-R/Microbiome2015_2020")

install.packages("tidyverse")

#Load packages
library(readxl)
library(ggplot2)
library(ggpubr)
library(rstatix);packageVersion("rstatix")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
library(tidyverse);packageVersion("tidyverse")
library("vegan"); packageVersion("vegan")

#Import Sample Data - "metadata"
meta <- read.csv("./data/Metadata_2015_2020_2.csv", h=T, sep = ",")

# Checking data structure
summary(meta)
str(meta)

#Edit order
all.data <- filter(meta, Dataset == "2020")
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("surface", "bottom"))
str(all.data)

#Edit colour pallete
location.col <- c("R-Estuary-1" = "#D55E00", 
                  "R-Estuary-2" = "#E69F00", 
                  "NS-Marine" = "#F0E442", 
                  "OS-Marine" = "#56B4E9",
                  "SM-Outfall" = "#009E73")
season.col <- c("winter" = "#56B4E9", 
                "spring" = "#F0E442", 
                "summer" = "#D55E00", 
                "autumn" = "#E69F00")


#Graph A: Temperature and dissolved oxygen
all.data.a <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = Temperature_sea, fill = Location), position = position_dodge2(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_point(mapping = aes(x = Season, y = Dissolved_Oxygen*4, shape = "Dissolved oxygen", group = Location), position = position_dodge2(width=.5), stat = 'identity') +
  scale_y_continuous(name = "Temperature [°C]", limits = c(0,35), 
                     sec.axis = sec_axis(~./4, name = expression (paste("DO [mg ", L^-1, "]")))) +
  scale_fill_manual(values=location.col)  +
  theme_bw() +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(vars(Depth)) +
  labs(shape = NULL, fill = NULL)
all.data.a 

#Graph B: Salinity
all.data.b <- ggplot(data = all.data) +
  geom_point(mapping = aes(x = Season, y = Salinity, colour = Location), position = position_dodge2(width=.7), stat = 'identity', show.legend = FALSE) +
  scale_y_continuous(name = "Salinity", limits = c(30,40)) +
  scale_colour_manual(values=location.col)   +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0,0,1), "lines"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(shape = NULL)
all.data.b 


#Graph C: DOC (y) / C/N (y2)
all.data.c <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = DOC, fill = Location), stat = 'identity', width = .5, position = position_dodge2(), show.legend = FALSE) +
  geom_point(mapping = aes(x = Season, y = C_N*9, shape = "C/N", group = Location), position = position_dodge2(width=.5), stat = 'identity') +
  scale_y_continuous(name = expression (paste("DOC [µmol", L^-1, "]")), limits = c(0,200),
                     sec.axis = sec_axis(~./9, name = expression (paste("C ", N^-1)))) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.c 

#Graph D: TDN
all.data.d <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = TDN, fill = Location), stat = 'identity', width = .5, position = position_dodge2(), show.legend = FALSE) +
  scale_y_continuous(name = expression (paste("TDN [µmol", L^-1, "]")), limits = c(0,30)) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "top", plot.margin = unit(c(0,0,0,1), "lines"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.d 


#Graph E: NH4 (y1) / NO2+NO3 (y2)
all.data.e <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = NH4, fill = Location), stat = 'identity', width = .5, position = position_dodge2(), show.legend = FALSE) +
  geom_point(aes(x = Season, y = NO2+NO3, shape = "NO2+NO3", group = Location), position = position_dodge2(width=.5), stat = 'identity') +
  scale_y_continuous(name = expression (paste("NH4 [µmol ", L^-1, "]")), limits = c(0,10),
                     sec.axis = sec_axis(~./1, name = expression (paste("NO2+NO3 [µmol ", L^-1, "]")))) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,1,0,0), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(shape = NULL, fill = NULL)
all.data.e


#Graph F: PO4 (y1) / SiO3 (y2)
all.data.f <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Season, y = PO4, fill = Location), stat = 'identity', width = .5, position = position_dodge2(), show.legend = FALSE) +
  geom_point(aes(x = Season, y = SiO3/5, shape = "SiO3", group = Location), position = position_dodge2(width=.5), stat = 'identity') +
  scale_y_continuous(name = expression (paste("PO4 [µmol ", L^-1, "]")), limits = c(0,3),
                     sec.axis = sec_axis(~.*5, name = expression (paste("SiO3 [µmol ", L^-1, "]")))) +
  scale_fill_manual(values=location.col) +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = c(0.9, 0.97), plot.margin = unit(c(0,0,0,1), "lines"), legend.background = element_rect(colour="black", linetype = "solid"), axis.text.x = element_blank(), axis.title.x = element_blank()) +
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
  theme(legend.position = "top", plot.margin = unit(c(0,1,0,0), "lines"), axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.g

#Graph H: Bacterial carbon production
all.data.h <- ggplot(data = all.data, aes(fill=Location)) +
  geom_bar(mapping = aes(x = Season, y = BCP), position = position_dodge(), stat = 'identity', width = .5, show.legend = FALSE) +
  geom_errorbar(mapping = aes(x = Season, ymin = BCP - BCP_STDEV, ymax = BCP + BCP_STDEV), position = position_dodge(width = .5), stat = 'identity', width = 0.25) +
  scale_y_continuous(name = expression(paste("BCP [N cell ", L^-1, " ", d^-1, "]"))) +
  scale_fill_manual(values=location.col)  +
  theme_bw() +
  facet_wrap(vars(Depth)) +
  theme(legend.position = "right", plot.margin = unit(c(0,0,0,1), "lines"), axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.h

ggarrange(all.data.a, all.data.b, all.data.c, all.data.d, all.data.e, all.data.f, all.data.g, all.data.h,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"), nrow = 4, ncol = 2, align = c("v"))


# CLEAN UP #################################################


# Clear plots
dev.off()  # But only if there IS a plot

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")  # ctrl+L

