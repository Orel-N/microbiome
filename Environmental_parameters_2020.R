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
meta <- read.csv("./data/Metadata_2015_2020.csv", h=T, sep = ",")

#Edit metadata
# Checking data structure
summary(meta)
str(meta)

meta[,'DIN']=round(meta[,'DIN'],1)
meta[,'Depth_m']=round(meta[,'Depth_m'],0)
meta[,'Coliform_bacteria']=round(meta[,'Coliform_bacteria'],0)
meta[,'Cyanobacteria_STDEV']=round(meta[,'Cyanobacteria_STDEV'],0)
meta[,'Cyanobacteria_STDEV']=round(meta[,'Cyanobacteria_STDEV'],0)
meta <- meta %>% mutate_at(12:18, round, 2)

meta$Dataset <- ifelse(meta$SampleID < 12, "2015", "2020")

#Save corrected table
write.csv(meta, "./data/Metadata_2015_2020.csv")


#Edit order
all.data <- filter(meta, Dataset == "2020")
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("surface", "bottom"))

all.data.a <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Temperature_sea, fill = "Temperature"), stat = 'identity', width = .5) +
  geom_point(mapping = aes(x = Location, y = Dissolved_Oxygen*4, shape = "Dissolved oxygen"), stat = 'identity') +
  scale_y_continuous(name = "Temperature [°C]", limits = c(0,35), 
                     sec.axis = sec_axis(~./4, name = expression (paste("DO [mg ", L^-1, "]")))) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.a 

ggsave("./output_graphs/a.pdf")

#Graph B: Salinity
all.data.b <- ggplot(data = all.data) +
  geom_point(mapping = aes(x = Location, y = Salinity, shape = "Salinity"), stat = 'identity') +
  scale_y_continuous(name = "Salinity", limits = c(30,45)) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL)
all.data.b 
ggsave("./output_graphs/b.pdf")

#Graph C: DOC/TDN (y1) C/N (y2)

#In order to bar plot side by side two variables (DOC,TDN), we have to do folowing data manipulation:
df <- all.data
c.df2 <- tidyr::pivot_longer(df, cols = c('DOC', 'TDN'), names_to = 'variable', 
                             values_to = "value")


all.data.c <- ggplot(data = c.df2) +
  geom_bar(mapping = aes(x = Location, y = value, fill = variable), stat = 'identity', width = .5, position = 'dodge') +
  geom_point(mapping = aes(x = Location, y = C_N*9, shape = "C/N"), stat = 'identity') +
  scale_y_continuous(name = expression (paste("DOC or TDN [µmol", L^-1, "]")), limits = c(0,200),
                     sec.axis = sec_axis(~./9, name = expression (paste("C ", N^-1)))) +
  scale_fill_manual(values = c("honeydew4", "honeydew3")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.c 
ggsave("./output_graphs/c.pdf")

#Graph D: DIN
d.df2 <- tidyr::pivot_longer(df, cols = c('NO2', 'NO3', 'NH4'), names_to = 'variable', 
                             values_to = "value")


all.data.d <- ggplot(data = d.df2) +
  geom_bar(mapping = aes(x = Location, y = value, fill = variable), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("[µmol ", L^-1, "]")), limits = c(0,10)) +
  scale_fill_manual(values = c("honeydew4", "honeydew3", "honeydew2")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.d
ggsave("./output_graphs/d.pdf")

#Graph E: PO4
all.data.e <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = PO43, fill = "PO43-"), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("[µmol ", L^-1, "]")), limits = c(0,1)) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.e
ggsave("./output_graphs/e.pdf")

#Graph E: SiO3
all.data.h <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = SiO3, fill = "SiO3"), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("[µmol ", L^-1, "]")), limits = c(0,20)) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.h
ggsave("./output_graphs/h.pdf")

#Graph F: Chl_A, Cyanobacteria
all.data.f <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Cyanobacteria, fill = "Cyanobacteria"), stat = 'identity', width = .5, position = 'dodge') +
  geom_errorbar(aes(x = Location, ymin = Cyanobacteria - Cyanobacteria_STDEV, ymax = Cyanobacteria + Cyanobacteria_STDEV), position = 'dodge', width = 0.25) +
  geom_point(mapping = aes(x = Location, y = Chl_A.1*100000000, shape = "Chl A"), stat = 'identity') +
  scale_y_continuous(name = expression (paste("Cyanobacteria [N cell", L^-1, "]")), limits = c(0,150000000),
                     sec.axis = sec_axis(~./100000000, name = expression (paste("Chl A [µg ", L^-1, "]")))) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.f
ggsave("./output_graphs/f.pdf")

#Graph G: Bacterial carbon production / bacterial abundance
all.data.g <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = BA_FC, fill = "Bacterial abundance"), stat = 'identity', width = .5, position = 'dodge') +
  geom_errorbar(aes(x = Location, ymin = BA_FC - BA_FC_STDEV, ymax = BA_FC + BA_FC_STDEV), position = 'dodge', width = 0.25) +
  geom_point(mapping = aes(x = Location, y = BCP*2, shape = "BCP"), stat = 'identity') +
  geom_errorbar(aes(x = Location, y = BCP*2, ymin = (BCP - BCP_STDEV)*2, ymax = (BCP + BCP_STDEV)*2), width = 0.25) +
  scale_y_continuous(name = expression(paste("Bacterial abundance [N cell ", L^-1, "]")),
                     sec.axis = sec_axis(~./2, name = expression (paste("BCP [N cell ", L^-1, " ", d^-1, "]")))) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(Depth~Season) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.g
ggsave("./output_graphs/g.pdf")

ggarrange(all.data.a, all.data.b, all.data.c, all.data.d, all.data.e, all.data.f,
          labels = c("A", "B", "C", "D", "E", "F"), nrow = 3, ncol = 2)

ggsave("./output_graphs/all.pdf")


# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")  # ctrl+L

