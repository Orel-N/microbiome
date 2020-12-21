#MICROBIOME - RDA ordination of bacterial community composition, dataset 2020
#Neža Orel, neza.orel@nib.si
#Script for correlation between env. parameters and ordination of bacterial community composition
############################################################################

#Set working directory
setwd("C:/Users/nezao/Documents/5-R/Microbiome2015_2020")

#Load packages
library("ggplot2"); packageVersion("ggplot2")
library(rstatix);packageVersion("rstatix")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
library(tidyverse);packageVersion("tidyverse")
library("vegan"); packageVersion("vegan")
library(RColorBrewer)
library(phyloseq)
library(ggpubr)

#Import collor palette
source("./scripts/Color_palettes.R")


####################################
#Import data
####################################

#Import Sample Data - "metadata"
meta <- read.csv("./data/Metadata_2015_2020.csv", h=T, sep = ",")

#Filter only 2020 dataset and edit table
all.data <- filter(meta, Dataset == "2020")
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("surface", "bottom"))
all.data$DistanceKp <- ifelse(all.data$Location == "R-Estuary-1", "750", ifelse(all.data$Location == "R-Estuary-2", "1500", 
                                                                                ifelse(all.data$Location == "NS-Marine", "2500", 
                                                                                       ifelse(all.data$Location == "OS-Marine", "16000", "18000"))))
all.data$DistanceKp <- as.numeric(all.data$DistanceKp)
all.data$DistanceWW <- ifelse(all.data$Location == "R-Estuary-1", "750", ifelse(all.data$Location == "R-Estuary-2", "1500", 
                                                                                ifelse(all.data$Location == "NS-Marine", "2500", 
                                                                                       ifelse(all.data$Location == "OS-Marine", "2600", "0"))))
all.data$DistanceWW <- as.numeric(all.data$DistanceWW)


#Import phyloseq (variance stabilased data)
phy_obj3 <- readRDS("./data/ps.vst.20.RDS")

#Import phyloseq (total abundance and relative abundance)
phy_obj <- readRDS("./data/phyloseqFiltered.RDS")
phy_obj <- subset_samples(phy_obj, Dataset == 2020)
phy_obj <- prune_taxa(taxa_sums(phy_obj)>0, phy_obj)
phy_obj.ra = transform_sample_counts(phy_obj, function (otu) {otu/sum(otu)})

#Edit phyloseq metadata: Add NO2+NO3
sample_data(phy_obj3)$NO2_NO3 <- sample_data(phy_obj3)$NO2 + sample_data(phy_obj3)$NO3
sample_data(phy_obj)$NO2_NO3 <- sample_data(phy_obj)$NO2 + sample_data(phy_obj)$NO3
sample_data(phy_obj3)$DistanceKp <- ifelse(sample_data(phy_obj3)$Location == "R-Estuary-1", "750", ifelse(sample_data(phy_obj3)$Location == "R-Estuary-2", "1500", 
                                                                                ifelse(sample_data(phy_obj3)$Location == "NS-Marine", "2500", 
                                                                                       ifelse(sample_data(phy_obj3)$Location == "OS-Marine", "16000", "18000"))))
sample_data(phy_obj3)$DistanceKp <- as.numeric(sample_data(phy_obj3)$DistanceKp)
sample_data(phy_obj3)$DistanceWW <- ifelse(sample_data(phy_obj3)$Location == "R-Estuary-1", "750", ifelse(sample_data(phy_obj3)$Location == "R-Estuary-2", "1500", 
                                                                                                          ifelse(sample_data(phy_obj3)$Location == "NS-Marine", "2500", 
                                                                                                                 ifelse(sample_data(phy_obj3)$Location == "OS-Marine", "2600", "0"))))
sample_data(phy_obj3)$DistanceWW <- as.numeric(sample_data(phy_obj3)$DistanceWW)

#Microbial indicators (from variance stabilased data)
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps.ind.vst <- subset_taxa(phy_obj3, Family %in% c(Mic_Ind$Family))
ps.ind.vst

#Most abundant taxa
most_abundant_taxa <- sort(taxa_sums(phy_obj3), TRUE) [1:2500]
ex2 <- prune_taxa(names(most_abundant_taxa), phy_obj3)

#Explore most abundant bacterial indicators
ps.ind <- subset_taxa(phy_obj.ra, Family %in% c(Mic_Ind$Family))
ind <- psmelt(ps.ind)
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ind,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]
top<-aggregate(Abundance~Family, ps_Mic_Ind_ra.melt.agg.genus, FUN=sum)


###############################################################
#Correlation between env. parameters - BEFORE FORWARD SELECTION
###############################################################

#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import env.par. data
metadata.raw <- all.data
env.par <- c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "TDN", "NO2_NO3", "PO4", "NH4", "C_N", "BA_FC", "BCP", "Chl_A", "SiO3", "DistanceKp", "DistanceWW")

#scale all parameters
metadata.scaled.SRF <- metadata.raw %>% 
  mutate_at(env.par, scale_par)%>%
  as.data.frame()

#correlation between the env parameters
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


###########################################
#RDA ordination - BEFORE FORWARD SELECTION
###########################################

#subset by fraction and remove NAs in metadata
phy_obj3
BAC_20.no.na <- phy_obj3 %>% 
  subset_samples(
    !is.na(Temperature_sea) &
      !is.na(Salinity) &
      !is.na(Dissolved_Oxygen_perc) &
      !is.na(DOC) &
      !is.na(TDN) & 
      !is.na(Chl_A.1) &
      !is.na(NO2) &
      !is.na(NO3) &      
      !is.na(PO4) &
      !is.na(NH4) &
      !is.na(SiO3) &
      Dataset == "2020")
BAC_20.no.na

##remove unobserved OTU and OTU from 2015 dataset
#BAC_20.no.na <- prune_taxa(taxa_sums(BAC_20.no.na)!=0,BAC_20.no.na)
#BAC_20.no.na

#extract and scale the env. parameters
BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "TDN", "PO4", "NO2_NO3", "NH4", "C_N", "DistanceWW", "Chl_A.1")]  
BAC_20.env <- as.data.frame(scale(BAC_20.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_20.otu <- t(otu_table(BAC_20.no.na))

#RDA analysis
BAC_20.rda.all <- rda (BAC_20.otu ~ ., data = BAC_20.env) # model including all variables


###########################################
#Forward selection of explanatory variables
###########################################

#based on the tutorial from David Zeleny Lab
#http://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel

BAC_20.rda.0 <- rda (BAC_20.otu ~ 1, data = BAC_20.env) # model containing only species matrix and intercept
BAC_20.rda.sel.os <- ordistep (BAC_20.rda.0, scope = formula (BAC_20.rda.all), direction = 'both') #stepwise selection

#Summary
BAC_20.rda.sel.os$anova

#variance partitioning 
varpart(BAC_20.otu, ~ Temperature_sea, ~ DOC, ~Dissolved_Oxygen, ~NH4, data = BAC_20.env[,c("Temperature_sea", "DOC", "Dissolved_Oxygen", "NH4")])


###############################################################
#Correlation between env. parameters - ONLY FORWARD SELECTED
###############################################################

#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import env.par. data
metadata.raw <- all.data
#env.par <- c("Temperature_sea", "DOC", "Dissolved_Oxygen_perc", "NH4")
#env.par <- c("Temperature_sea","Salinity", "DOC", "NO2_NO3")
env.par <- c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "NO2_NO3", "PO4", "NH4", "DistanceKp")

#scale all parameters
metadata.scaled.SRF <- metadata.raw %>% 
  mutate_at(env.par, scale_par)%>%
  as.data.frame()

#correlation between the env parameters
envpar_corr <- metadata.scaled.SRF %>% 
  select(env.par)%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr.pvalues<- envpar_corr %>% cor_get_pval()

#plot
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)


##########################################################################
#RDA ordination of the fitted model of bacterial composition data 
#constrained by environmental variables - ONLY FORWARD SELECTED ENV. VAR
########################################################################

#subset by fraction and remove NAs in metadata
phy_obj3
BAC_20.no.na <- phy_obj3 %>% 
  subset_samples(
    !is.na(Temperature_sea) &
      !is.na(Salinity) &
      !is.na(Dissolved_Oxygen) &
      !is.na(DOC) &
      #!is.na(TDN) & 
      #!is.na(Chl_A.1) &
      !is.na(NO2_NO3) &
      #!is.na(NO3) &      
      !is.na(PO4) &
      !is.na(NH4) &
      #!is.na(SiO3) &
      Dataset == "2020")
BAC_20.no.na

##remove unobserved OTU and OTU from 2015 dataset
#BAC_20.no.na <- prune_taxa(taxa_sums(BAC_20.no.na)!=0,BAC_20.no.na)
#BAC_20.no.na

#extract and scale the env. parameters - only forward selected
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea", "DOC", "Dissolved_Oxygen", "NH4")]  
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity", "DOC", "NO2_NO3")]  
BAC_20.env <- data.frame(sample_data(BAC_20.no.na)) [c("Temperature_sea","Dissolved_Oxygen", "DOC", "PO4", "NH4", "DistanceKp")]
BAC_20.env <- as.data.frame(scale(BAC_20.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_20.otu <- t(otu_table(BAC_20.no.na))

#RDA analysis
BAC_20.rda.all <- rda (BAC_20.otu ~ ., data = BAC_20.env) # model including all variables

#The summary of RDA
summary(BAC_20.rda.all, display = NULL)
anova.cca(BAC_20.rda.all, step=1000)
#anova.cca(BAC_20.rda.all, step=1000, by="axis")

#generate an RDA plot ##########################################

#Extract RDA scores
BAC_20.rda.scores <- vegan::scores(BAC_20.rda.all,display=c("sp","wa","lc","bp","cn"))

#Make data frame with scores - Locations
BAC_20.rda.sites <- data.frame(BAC_20.rda.scores$sites)
BAC_20.rda.sites$SampleID <- as.character(rownames(BAC_20.rda.sites))
sample_data(BAC_20.no.na)$SampleID <- as.character(rownames(sample_data(BAC_20.no.na)))
BAC_20.rda.sites <- BAC_20.rda.sites %>%
  left_join(sample_data(BAC_20.no.na))

#Make data frame with scores - Species
BAC_20.rda.species <- data.frame(BAC_20.rda.scores$species)
BAC_20.rda.species$species_names <- rownames(BAC_20.rda.species)
TAX <- data.frame(as(tax_table(BAC_20.no.na), "matrix"))
TAX$species_names <- rownames(TAX)
BAC_20.rda.species <- BAC_20.rda.species %>%
  left_join(TAX)

#Select different groups for visualization
BAC_20.rda.species <- filter(BAC_20.rda.species, species_names %in% rownames(otu_table(ex2)))
MIC_20.rda.species <- filter(BAC_20.rda.species, Family %in% c(Mic_Ind$Family))
Aero.rda.species <- filter(BAC_20.rda.species, Family == "Aeromonadaceae")
Arcobacter.rda.species <- filter(BAC_20.rda.species, Family == "Arcobacteraceae")
Bacter.rda.species <- filter(BAC_20.rda.species, Family == "Bacteroidaceae")
Enterbac.rda.species <- filter(BAC_20.rda.species, Family == "Enterobacteriaceae")
Enter.rda.species <- filter(BAC_20.rda.species, Family == "Enterococcaceae")
Lach.rda.species <- filter(BAC_20.rda.species, Family == "Lachnospiraceae")
Moraxellaceae.rda.species <- filter(BAC_20.rda.species, Family == "Moraxellaceae")
Pseu.rda.species <- filter(BAC_20.rda.species, Family == "Pseudomonadaceae")
Rum.rda.species <- filter(BAC_20.rda.species, Family == "Ruminococcaceae")
Vibrio.rda.species <- filter(BAC_20.rda.species, Family == "Vibrionaceae")
Mycobact.rda.species <- filter(BAC_20.rda.species, Family == "Mycobacteriaceae")
Clostridia.rda.species <- filter(BAC_20.rda.species, Family == "Clostridiaceae")
Bdelovib.rda.species <- filter(BAC_20.rda.species, Family == "Bdellovibrionaceae")


c.1 <- filter(BAC_20.rda.species, Class == "Alphaproteobacteria")
c.2 <- filter(BAC_20.rda.species, Class == "Bacteroidia")
c.3 <- filter(BAC_20.rda.species, Class == "Campylobacteria")
c.4 <- filter(BAC_20.rda.species, Class == "Cyanobacteriia")
c.5 <- filter(BAC_20.rda.species, Class == "Gammaproteobacteria")

c.6 <- filter(BAC_20.rda.species, Order == "Rhodobacterales")
c.7 <- filter(BAC_20.rda.species, Order == "SAR11_clade")
c.8 <- filter(BAC_20.rda.species, Order == "Thiomicrospirales")
c.9 <- filter(BAC_20.rda.species, Order == "Cellvibrionales")
c.10 <- filter(BAC_20.rda.species, Order == "Alteromonadales")
c.11 <- filter(BAC_20.rda.species, Order == "Flavoacteriales")


#BAC_20.rda.species <- BAC_20.rda.species[abs(BAC_20.rda.species$RDA1)>0.5 | abs(BAC_20.rda.species$RDA2)>0.5,]
#BAC_20.rda.species <- BAC_20.rda.species[abs(BAC_20.rda.species$RDA1)>0.35 | abs(BAC_20.rda.species$RDA2)>0.35,]


#Draw biplots
BAC_20.rda.arrows<- BAC_20.rda.scores$biplot*10
colnames(BAC_20.rda.arrows)<-c("x","y")
BAC_20.rda.arrows <- as.data.frame(BAC_20.rda.arrows)
BAC_20.rda.evals <- 100 * summary(BAC_20.rda.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot A: Sites and env.par
BAC_20.rda.plot <- ggplot() +
  geom_point(data = BAC_20.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Depth), 
            size = 4) +
  geom_text(data = BAC_20.rda.sites,aes(x = RDA1, y = RDA2,label = Location), 
            nudge_y= -0.3,size=2) +
  scale_colour_manual(values=season.col) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic()+
  theme(legend.position = "top")+
  geom_segment(data=BAC_20.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_20.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_20.rda.arrows)),color="black",alpha=0.5)
BAC_20.rda.plot

#Plot B: Species (phylum) and env.par
BAC_20.rda.B.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Phylum), #alpha = 1/10,
              size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme(legend.position = "top")+
  theme_classic()+
  geom_segment(data=BAC_20.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_20.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_20.rda.arrows)),color="black",alpha=0.5)
BAC_20.rda.B.plot

#Plot C: Mic_ind
Mic_ind.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=MIC_20.rda.species, aes(x = RDA1*15, y = RDA2*15, fill = "Microbial indicators"), color="red", #show.legend = FALSE,
          size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank())
#+
  #scale_colour_manual(values=indicators.col)
Mic_ind.plot


#Plot D: Separeted families
m.1.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Arcobacter.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.1.plot

m.2.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Aero.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.2.plot

m.3.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Bacter.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.3.plot

m.4.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Enterbac.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.4.plot

m.13.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Enter.rda.species , aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.13.plot

m.5.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Lach.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic()+
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.5.plot

m.6.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Moraxellaceae.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.6.plot

m.7.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Pseu.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.7.plot

m.8.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Rum.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.8.plot

m.9.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Vibrio.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.9.plot

m.10.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Mycobact.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col)
m.10.plot

m.11.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Clostridia.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.11.plot

m.12.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=Bdelovib.rda.species, aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)
m.12.plot

ggarrange(BAC_20.rda.plot, Mic_ind.plot, m.9.plot, m.1.plot, m.12.plot, m.4.plot, ncol=2, nrow=3)

###############################################
#Plot D:
c.1.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.1, aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col)
c.1.plot

c.2.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.2, aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col)
c.2.plot

c.3.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.3, aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col)
c.3.plot

c.4.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.4, aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col)
c.4.plot

c.5.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.5, aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=phyla.col)
c.5.plot

ggarrange(BAC_20.rda.plot, c.1.plot, c.5.plot, c.2.plot, c.4.plot, c.3.plot, ncol=2, nrow=3)

###
c.6.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.6, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.6.plot

c.7.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.7, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.7.plot

c.8.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.8, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.8.plot

c.9.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.9, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.9.plot

c.10.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.10, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.10.plot

c.11.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=c.11, aes(x = RDA1*15, y = RDA2*15, color = Order),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values="red")
c.11.plot

# CLEAN UP #################################################

# Clear packages
detach("package:phyloseq", unload=TRUE)
detach("package:rstatix", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")  # ctrl+L