#MICROBIOME - RDA ordination of bacterial community composition
#Neža Orel, neza.orel@nib.si
#Script for correlation between env. parameters and RDA ordination of bacterial community composition
#Data was collected during in situ survey 2020, project MIKROBIOM

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

#Import Sample Data - "metadata"
meta <- read.csv("./data/Metadata_2015_2020.csv", h=T, sep = ",")

#Filter only 2020 dataset and edit table
all.data <- filter(meta, Dataset == "2020")
all.data$Season <- factor(all.data$Season, levels = c("winter", "spring", "summer", "autumn"))
all.data$Location <- factor(all.data$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
all.data$Depth <- factor(all.data$Depth, levels = c("surface", "bottom"))

#Import phyloseq (variance stabilased data)
phy_obj3 <- readRDS("./data/ps.vst.RDS")
phy_obj3

#Microbial indicators (from variance stabilased data)
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps.ind.vst <- subset_taxa(phy_obj3, Family %in% c(Mic_Ind$Family))
ps.ind.vst

################################
#Correlation between env. parameters
################################

#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import env.par. data
metadata.raw <- all.data
#env.par <- c("Temperature_sea","Salinity","Chl_A.1", "Dissolved_Oxygen_perc", "DOC", "TDN", "NO2", "NO3", "PO4", "SiO3", "NH4")
env.par <- c("Temperature_sea","Salinity","DOC", "Dissolved_Oxygen_perc", "DIN")
#env.par <- c("Temperature_sea","Salinity","DOC", "DIN")

#scale all the parameters
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

#####################################
#RDA ordination of the fitted model of bacterial composition data constrained by environmental variables
####################################

#subset by fraction and remove NAs in metadata
BAC_20.no.na <- phy_obj3 %>% 
  subset_samples(
    !is.na(Temperature_sea) &
      !is.na(Salinity) &
      !is.na(Dissolved_Oxygen_perc) &
      !is.na(DOC) &
      !is.na(TDN) & 
      #      !is.na(Chl_A.1) &
      !is.na(NO2) &
      !is.na(NO3) &      
      !is.na(PO4) &
      !is.na(NH4) &
      !is.na(SiO3) &
      Dataset == "2020")
BAC_20.no.na

##remove unobserved OTU and OTU from 2015 dataset
BAC_20.no.na <- prune_taxa(taxa_sums(BAC_20.no.na)>0,BAC_20.no.na)
BAC_20.no.na

#extract and scale the env. parameters
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity","Chl_A.1", "Dissolved_Oxygen_perc", "DOC", "TDN", "NO2", "NO3", "PO4", "SiO3", "NH4")]  
BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity","DOC", "Dissolved_Oxygen_perc", "DIN")]  
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity","DOC", "DIN")]  

BAC_20.env <- as.data.frame(scale(BAC_20.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_20.otu <- t(otu_table(BAC_20.no.na))

#RDA analysis
BAC_20.rda.all <- rda (BAC_20.otu ~ ., data = BAC_20.env) # model including all variables 

#generate an RDA plot 
BAC_20.rda.scores <- vegan::scores(BAC_20.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC_20.rda.sites <- data.frame(BAC_20.rda.scores$sites)
BAC_20.rda.sites$SampleID <- as.character(rownames(BAC_20.rda.sites))
sample_data(BAC_20.no.na)$SampleID <- as.character(rownames(sample_data(BAC_20.no.na)))
BAC_20.rda.sites <- BAC_20.rda.sites %>%
  left_join(sample_data(BAC_20.no.na))

#Draw biplots
BAC_20.rda.arrows<- BAC_20.rda.scores$biplot*5
colnames(BAC_20.rda.arrows)<-c("x","y")
BAC_20.rda.arrows <- as.data.frame(BAC_20.rda.arrows)
BAC_20.rda.evals <- 100 * (BAC_20.rda.all$CCA$eig / sum(BAC_20.rda.all$CCA$eig))

#Plot 
BAC_20.rda.plot <- ggplot() +
  geom_point(data = BAC_20.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Depth), 
             size = 4) +
  geom_text(data = BAC_20.rda.sites,aes(x = RDA1, y = RDA2,label = Location), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  #scale_x_reverse()+ 
  theme(legend.position = "top")+
  geom_segment(data=BAC_20.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_20.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_20.rda.arrows)),color="black",alpha=0.5)
BAC_20.rda.plot


#####################################
#RDA ordination of the fitted model of bacterial composition data constrained by environmental variables - Microbial indicators
####################################

#subset by fraction and remove NAs in metadata
Mic_20.no.na <- ps.ind.vst %>% 
  subset_samples(
    !is.na(Temperature_sea) &
      !is.na(Salinity) &
      !is.na(Dissolved_Oxygen_perc) &
      !is.na(DOC) &
      !is.na(TDN) & 
      #      !is.na(Chl_A.1) &
      !is.na(NO2) &
      !is.na(NO3) &      
      !is.na(PO4) &
      !is.na(NH4) &
      !is.na(SiO3) &
      Dataset == "2020")
Mic_20.no.na

##remove unobserved OTU and OTU from 2015 dataset
Mic_20.no.na <- prune_taxa(taxa_sums(Mic_20.no.na)>0, Mic_20.no.na)
Mic_20.no.na

#extract and scale the env. parameters
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity","Chl_A.1", "Dissolved_Oxygen_perc", "DOC", "TDN", "NO2", "NO3", "PO4", "SiO3", "NH4")]  
Mic_20.env <- data.frame(sample_data(Mic_20.no.na))[c("Temperature_sea","Salinity","DOC", "Dissolved_Oxygen_perc", "DIN")]  
#BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity","DOC", "DIN")]  

Mic_20.env <- as.data.frame(scale(Mic_20.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
Mic_20.otu <- t(otu_table(Mic_20.no.na))

#RDA analysis
Mic_20.rda.all <- rda (Mic_20.otu ~ ., data = Mic_20.env) # model including all variables 

#generate an RDA plot 
Mic_20.rda.scores <- vegan::scores(Mic_20.rda.all,display=c("sp","wa","lc","bp","cn"))
Mic_20.rda.sites <- data.frame(Mic_20.rda.scores$sites)
Mic_20.rda.sites$SampleID <- as.character(rownames(Mic_20.rda.sites))
sample_data(Mic_20.no.na)$SampleID <- as.character(rownames(sample_data(Mic_20.no.na)))
Mic_20.rda.sites <- Mic_20.rda.sites %>%
  left_join(sample_data(Mic_20.no.na))

#Draw biplots
Mic_20.rda.arrows<- Mic_20.rda.scores$biplot*10
colnames(Mic_20.rda.arrows)<-c("x","y")
Mic_20.rda.arrows <- as.data.frame(Mic_20.rda.arrows)
Mic_20.rda.evals <- 100 * (Mic_20.rda.all$CCA$eig / sum(Mic_20.rda.all$CCA$eig))

#Plot 
Mic_20.rda.plot <- ggplot() +
  geom_point(data = Mic_20.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Depth), 
             size = 4) +
  geom_text(data = Mic_20.rda.sites,aes(x = RDA1, y = RDA2,label = Location), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(Mic_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Mic_20.rda.evals[2], 2))) +
  #scale_x_reverse()+ 
  theme(legend.position = "top")+
  geom_segment(data=Mic_20.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(Mic_20.rda.arrows*1.2),
            aes(x, y, label = rownames(Mic_20.rda.arrows)),color="black",alpha=0.5)
Mic_20.rda.plot



#####################################
#Forward selection of explanatory variables
#####################################
#based on the tutorial from David Zeleny Lab
#http://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel

BAC_20.rda.0 <- rda (BAC_20.otu ~ 1, data = BAC_20.env) # model containing only species matrix and intercept
BAC_20.rda.sel.os <- ordistep (BAC_20.rda.0, scope = formula (BAC_20.rda.all), direction = 'both') #stepwise selection

BAC_20.rda.sel.os$anova

#variance partitioning 
varpart(BAC_20.otu, ~ Salinity, ~ Temperature_sea, ~ DOC, data = BAC_20.env[,c("Salinity", "Temperature_sea", "DOC")])

# CLEAN UP #################################################

# Clear packages
detach("package:phyloseq", unload=TRUE)
detach("package:rstatix", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")  # ctrl+L
