###################################################
#RDA ordination of bacterial community composition
###################################################

#Load packages
library(ggplot2); packageVersion("ggplot2")
library(dplyr);packageVersion("dplyr")
library(phyloseq)
library(tidyr);packageVersion("tidyr")
library(rstatix);packageVersion("rstatix")
library(vegan); packageVersion("vegan")
library(ggpubr)

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


############
#Import data
############

#Import Sample Data - "metadata", filter only 2020 dataset and edit table
metadata.raw <- read.csv("./data/Metadata_2020.csv", sep = ",")

metadata.raw$Season <- factor(metadata.raw$Season, levels = c("winter", "spring", "summer", "autumn"))
metadata.raw$Location <- factor(metadata.raw$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
metadata.raw$Depth <- factor(metadata.raw$Depth, levels = c("surface", "bottom"))

#Import phyloseq (variance stabilased data)
phy_obj3 <- readRDS("./data/phyloseqVarStab20.RDS")

#Import phyloseq (total abundance), calculate relative abundance and subset 2020 dataset
phy_obj <- readRDS("./data/phyloseqPrevFiltered20.RDS")

phy_obj.ra = transform_sample_counts(phy_obj, function (otu) {otu/sum(otu)})

#Microbial indicators (from variance stabilased data)
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps.ind.vst <- subset_taxa(phy_obj3, Family %in% c(Mic_Ind$Family))
ps.ind.vst


#######################################################
#Correlation between env. parameters - BEFORE SELECTION
#######################################################

#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import env.par. data
env.par <- c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "TDN", "NO2_NO3", "PO4", "NH4", "C_N", "BA_FC", "BCP_C", "Chl_A.1", "SiO3")

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


##################################
#RDA ordination - BEFORE SELECTION
##################################

#extract and scale the env. parameters
BAC_20.env <- data.frame(sample_data(BAC_20.no.na))[c("Temperature_sea","Salinity", "Dissolved_Oxygen", "DOC", "TDN", "PO4", "NO2_NO3", "NH4", "C_N")]  
BAC_20.env <- as.data.frame(scale(BAC_20.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_20.otu <- t(otu_table(BAC_20.no.na))

#RDA analysis
BAC_20.rda.all <- rda (BAC_20.otu ~ ., data = BAC_20.env) # model including all variables


###########################################
#Selection of explanatory variables
###########################################

#based on the tutorial from David Zeleny Lab
#http://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel

BAC_20.rda.0 <- rda (BAC_20.otu ~ 1, data = BAC_20.env) # model containing only species matrix and intercept
BAC_20.rda.sel.os <- ordistep (BAC_20.rda.0, scope = formula (BAC_20.rda.all), direction = 'both',
                               permutations = how(nperm = 999), steps = 100) #stepwise selection

#summary of RDA
summary(BAC_20.rda.sel.os, display = NULL)
vif.cca(BAC_20.rda.sel.os)

RsquareAdj(BAC_20.rda.sel.os)
anova(BAC_20.rda.sel.os)
anova(BAC_20.rda.sel.os,by="terms")
anova(BAC_20.rda.sel.os,by="axis")

#Variance partitioning 
varpart(BAC_20.otu, ~Temperature_sea, ~Dissolved_Oxygen, ~DOC, ~PO4,  data = BAC_20.env[,c("Temperature_sea", "Dissolved_Oxygen", "DOC", "PO4")])


####################
#Figure 5: RDA plot 
####################

#Extract RDA scores
BAC_20.rda.scores <- vegan::scores(BAC_20.rda.sel.os,display=c("sp","wa","lc","bp","cn"))

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

MIC_20.rda.species <- filter(BAC_20.rda.species, Family %in% c(Mic_Ind$Family))

#Biplot arrows
BAC_20.rda.arrows<- BAC_20.rda.scores$biplot*8
colnames(BAC_20.rda.arrows)<-c("x","y")
BAC_20.rda.arrows <- as.data.frame(BAC_20.rda.arrows)
BAC_20.rda.evals <- 100 * summary(BAC_20.rda.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot A: Sites and env.par
BAC_20.rda.plot <- ggplot() +
  geom_point(data = BAC_20.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Depth), 
            size = 3) +
  scale_colour_manual(values=season.col) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic()+
  theme(legend.position = "top", legend.title = element_text(size=10), legend.text = element_text(size=10))+
  geom_segment(data=BAC_20.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_20.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_20.rda.arrows)),color="black",alpha=1, size = 4)+
  guides(color = guide_legend(override.aes = list(size = 3)) )

BAC_20.rda.plot2 <- BAC_20.rda.plot + scale_shape(guide=FALSE) + coord_fixed() + xlim(-10, 10)
BAC_20.rda.plot2

ggsave("./output_graphs/RDA_samples.pdf", BAC_20.rda.plot)

#Plot B:
c.1.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Class == "Alphaproteobacteria"), aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=phyla.col) + coord_fixed() 
c.1.plot

#Plot C:
c.2.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Class == "Gammaproteobacteria"), aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col) + coord_fixed() 
c.2.plot

#Plot D:
c.3.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Class == "Bacteroidia"), aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col) + coord_fixed()
c.3.plot

#Plot E:
c.4.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Class == "Cyanobacteriia"), aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=phyla.col) + coord_fixed()
c.4.plot

#Plot F:
c.5.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Class == "Campylobacteria"), aes(x = RDA1*15, y = RDA2*15, color = Class),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=phyla.col) + coord_fixed() 
c.5.plot

#Plot G: Mic_ind
Mic_ind.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=MIC_20.rda.species, aes(x = RDA1*15, y = RDA2*15, fill = "Microbial indicators"), color="orange", #show.legend = FALSE,
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank()) + coord_fixed() 
Mic_ind.plot

ggarrange(c.1.plot, c.2.plot, c.3.plot, c.4.plot, c.5.plot, Mic_ind.plot, ncol=2, nrow=3, labels = c("B", "C", "D", "E", "F", "G"))
ggsave("./output_graphs/RDA_DominantClass.pdf", last_plot())


##########################################
#Supplementary figure 4: RDA plot Mic. Ind
##########################################

#Plot C2: Enterobacteriaceae
m.1.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Family == "Enterobacteriaceae"), aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col) + coord_fixed() 
m.1.plot

#Plot C2: Bacteroidaceae
m.2.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Family == "Bacteroidaceae"), aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col) + coord_fixed()
m.2.plot

#Plot C2: Lachnospiraceae
m.3.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Family == "Lachnospiraceae"), aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col) + coord_fixed()
m.3.plot

#Plot D2: Vibrionaceae
m.4.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Family == "Vibrionaceae"), aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_colour_manual(values=indicators.col)+ coord_fixed()
m.4.plot

#Plot E2: Arcobacteraceae
m.5.plot <- ggplot() +
  geom_point(data=BAC_20.rda.species, aes(x = RDA1*15, y = RDA2*15), alpha = 1/15,
             size=2) +
  geom_point(data=subset(BAC_20.rda.species, Family == "Arcobacteraceae"), aes(x = RDA1*15, y = RDA2*15, color = Family),
             size=2) +
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_20.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_20.rda.evals[2], 2))) +
  theme_classic() +
  theme(legend.position = "top")+
  scale_colour_manual(values=indicators.col) + coord_fixed()
m.5.plot

ggarrange(Mic_ind.plot, m.1.plot, m.2.plot, m.3.plot, m.4.plot, m.5.plot, labels = c("B", "C", "D", "E", "F", "G"), ncol=2, nrow=3)

ggsave("./output_graphs/RDA_MicInd.pdf", last_plot())


# CLEAN UP #################################################

# Clear packages
detach("package:phyloseq", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)

# Clear plots
dev.off()  

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
