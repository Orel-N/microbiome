#MICROBIOME - Dynamics of water biotic and abiotic paremeters
#Neža Orel, neza.orel@nib.si
#Script for graphs of accompining physical, chemical and biological properties.
#Data was collected during in situ survey 2015, project MIKROBIOM

############################################################################

#Set working directory
setwd("C:/Users/nezao/Documents/5-R/Microbiome2015_2020")

#Load packages
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)

#Import data
all.data <- read_excel("data/EnvPar_2015.xlsx")
View(all.data)

# Checking data structure
summary(all.data)
str(all.data)

#Edit order
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Mouth", "R-Estuary-1", "NS-Marine"))


all.data.a <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = Temperature_sea, fill = "Temperature"), stat = 'identity', width = .5) +
  geom_point(mapping = aes(x = Location, y = Dissolved_Oxygen*3.2, shape = "Dissolved oxygen"), stat = 'identity') +
  scale_y_continuous(name = "Temperature [°C]", limits = c(0,35), 
                     sec.axis = sec_axis(~./3.2, name = expression (paste("DO [mg ", L^-1, "]")))) +
  scale_fill_manual(values = c("honeydew4")) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_grid(~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.a 

ggsave("./output_graphs/a15.pdf")

#Graph B: Salinity
all.data.b <- ggplot(data = all.data) +
  geom_point(mapping = aes(x = Location, y = Salinity, shape = "Salinity"), stat = 'identity') +
  scale_y_continuous(name = "Salinity") +
  scale_fill_manual(values = c("honeydew4")) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_grid(~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL)
all.data.b 
ggsave("./output_graphs/b15.pdf")

#Graph C: DOC
all.data.c <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = DOC, fill = "DOC"), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("DOC [µmol", L^-1, "]")), limits = c(0,200),
                     sec.axis = sec_axis(~./9, name = expression (paste("C ", N^-1)))) +
  scale_fill_manual(values = c("honeydew4")) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_grid(~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.c 
ggsave("./output_graphs/c15.pdf")

#Graph D: DIN
df <- all.data
d.df2 <- tidyr::pivot_longer(df, cols = c('NO2', 'NO3', 'NH4'), names_to = 'variable', 
                             values_to = "value")

all.data.d <- ggplot(data = d.df2) +
  geom_bar(mapping = aes(x = Location, y = value, fill = variable), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("[µmol ", L^-1, "]"))) +
  scale_fill_manual(values = c("honeydew4", "honeydew3", "honeydew2")) +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  facet_grid(~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.d
ggsave("./output_graphs/d15.pdf")

#Graph E: PO4
all.data.e <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = PO43, fill = "PO43-"), stat = 'identity', width = .5, position = 'dodge') +
  scale_y_continuous(name = expression (paste("[PO43 µmol ", L^-1, "]")), limits = c(0,1)) +
  scale_fill_manual(values = c("honeydew4")) +
  theme_bw() +
  facet_grid(~Season, scales = "free") +
  scale_x_discrete(drop = FALSE) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.e
ggsave("./output_graphs/e15.pdf")

#Graph G: Bacterial carbon production / bacterial abundance
all.data.g <- ggplot(data = all.data) +
  geom_bar(mapping = aes(x = Location, y = BA, fill = "Bacterial abundance"), stat = 'identity', width = .5, position = 'dodge') +
  geom_errorbar(aes(x = Location, ymin = BA - BA_STDEV, ymax = BA + BA_STDEV), position = 'dodge', width = 0.25) +
  scale_y_continuous(name = expression(paste("Bacterial abundance [N cell ", L^-1, "]"))) +
  scale_fill_manual(values = c("honeydew4")) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_grid(~Season, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(shape = NULL, fill = NULL)
all.data.g
ggsave("./output_graphs/g15.pdf")

ggarrange(all.data.a, all.data.b, all.data.c, all.data.d, all.data.e, all.data.g, 
          labels = c("A", "B", "C", "D", "E", "G"), nrow = 3, ncol = 2)

ggsave("./output_graphs/all_15.pdf")

# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")  # ctrl+L

