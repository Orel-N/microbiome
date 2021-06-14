########################################
#Shared ASVs between different Locations
########################################

#Load packages
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(dplyr)
library(nVennR)
library(RColorBrewer)

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


############
#Import data
############

#Script for presence/absence matrix
source("./scripts/pres_abs_matrix.R")

#Import data
BAC_pruned <- readRDS("./data/phyloseqPrevFiltered20.RDS")
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")


##########################################################
#Diagram of bacterial OTU overlap between different groups 
##########################################################

#Subset by location
BAC_ERI2 <- subset_samples(BAC_pruned, Location == "R-Estuary-1")
BAC_ERI2 <- prune_taxa(taxa_sums(BAC_ERI2)>0, BAC_ERI2)
BAC_0014 <- subset_samples(BAC_pruned, Location == "R-Estuary-2")
BAC_0014 <- prune_taxa(taxa_sums(BAC_0014)>0, BAC_0014)
BAC_000K <- subset_samples(BAC_pruned, Location == "NS-Marine")
BAC_000K <- prune_taxa(taxa_sums(BAC_000K)>0, BAC_000K)
BAC_00BF <- subset_samples(BAC_pruned, Location == "OS-Marine")
BAC_00BF <- prune_taxa(taxa_sums(BAC_00BF)>0, BAC_00BF)
BAC_CN01 <- subset_samples(BAC_pruned, Location == "SM-Outfall")
BAC_CN01 <- prune_taxa(taxa_sums(BAC_CN01)>0, BAC_CN01)

#Mic indicators
Mic_pruned <- subset_taxa(BAC_pruned, Family %in% c(Mic_Ind$Family))
Mic_pruned <- prune_taxa(taxa_sums(Mic_pruned)>0, Mic_pruned)
Mic_ERI2 <- subset_samples(Mic_pruned, Location == "R-Estuary-1")
Mic_ERI2 <- prune_taxa(taxa_sums(Mic_ERI2)>0, Mic_ERI2)
Mic_0014 <- subset_samples(Mic_pruned, Location == "R-Estuary-2")
Mic_0014 <- prune_taxa(taxa_sums(Mic_0014)>0, Mic_0014)
Mic_000K <- subset_samples(Mic_pruned, Location == "NS-Marine")
Mic_000K <- prune_taxa(taxa_sums(Mic_000K)>0, Mic_000K)
Mic_00BF <- subset_samples(Mic_pruned, Location == "OS-Marine")
Mic_00BF <- prune_taxa(taxa_sums(Mic_00BF)>0, Mic_00BF)
Mic_CN01 <- subset_samples(Mic_pruned, Location == "SM-Outfall")
Mic_CN01 <- prune_taxa(taxa_sums(Mic_CN01)>0, Mic_CN01)

#make a list
y <- list()
y[["REstuary1"]] <- as.character(row.names(otu_table(BAC_ERI2)))
y[["REstuary2"]] <- as.character(row.names(otu_table(BAC_0014)))
y[["NSMarine"]] <- as.character(row.names(otu_table(BAC_000K)))
y[["OSMarine"]] <- as.character(row.names(otu_table(BAC_00BF)))
y[["SMOutfall"]] <- as.character(row.names(otu_table(BAC_CN01)))

y.mic <- list()
y.mic[["REstuary1"]] <- as.character(row.names(otu_table(Mic_ERI2)))
y.mic[["REstuary2"]] <- as.character(row.names(otu_table(Mic_0014)))
y.mic[["NSMarine"]] <- as.character(row.names(otu_table(Mic_000K)))
y.mic[["OSMarine"]] <- as.character(row.names(otu_table(Mic_00BF)))
y.mic[["SMOutfall"]] <- as.character(row.names(otu_table(Mic_CN01)))

#generate presence/abscence matrix
otu_overlaps <- pres_abs_matrix(y)    
otu_overlaps$OTU <- rownames(otu_overlaps)

otu_overlaps.mic <- pres_abs_matrix(y.mic)    
otu_overlaps.mic$OTU <- rownames(otu_overlaps.mic)


########################
#Explore overlaping OTU
########################

taxonomy <- as.data.frame(tax_table(BAC_pruned))
taxonomy$OTU <- rownames(taxonomy)

taxonomy.mic <- as.data.frame(tax_table(Mic_pruned))
taxonomy.mic$OTU <- rownames(taxonomy.mic)

otu_overlaps_merged <- full_join(taxonomy,otu_overlaps, by = c("OTU"))
otu_overlaps_merged.mic <- full_join(taxonomy.mic,otu_overlaps.mic, by = c("OTU"))

otu_overlaps_merged_B <- otu_overlaps_merged
otu_overlaps_merged.mic_B <- otu_overlaps_merged.mic


#####################################
#Relative abundance of overlaping OTU
#####################################

#add abundance in each fraction
#transform data
BAC_pruned.ra <- transform_sample_counts(BAC_pruned, function(x) x / sum(x))
Mic_pruned.ra <- subset_taxa(BAC_pruned.ra, Family %in% c(Mic_Ind$Family))
Mic_pruned.ra <- prune_taxa(taxa_sums(Mic_pruned.ra)>0, Mic_pruned.ra)

#calculate mean abundance for each ASV
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)
BAC_pruned.ra.long.agg <- aggregate(Abundance~OTU, BAC_pruned.ra.long, FUN = mean)
BAC_pruned.ra.long.agg$Abundance <- BAC_pruned.ra.long.agg$Abundance*100

Mic_pruned.ra.long <- psmelt(Mic_pruned.ra)
Mic_pruned.ra.long.agg <- aggregate(Abundance~OTU, Mic_pruned.ra.long, FUN = mean)
Mic_pruned.ra.long.agg$Abundance <- Mic_pruned.ra.long.agg$Abundance*100

#Join otu table of overlaping OTU and abundance table
otu_overlaps_merged <- full_join(otu_overlaps_merged, BAC_pruned.ra.long.agg, by = "OTU")
otu_overlaps_merged.mic <- full_join(otu_overlaps_merged.mic, Mic_pruned.ra.long.agg, by = "OTU")

#set metadata
sets <- names(y)
metadata <- as.data.frame(sets)

sets.mic <- names(y.mic)
metadata.mic <- as.data.frame(sets.mic)


#######################
#Figure 6: Venn diagram
#######################

#All
RE1 <- subset(otu_overlaps, REstuary1 == "1")$OTU
RE2 <- subset(otu_overlaps, REstuary2 == "1")$OTU
NSM <- subset(otu_overlaps, NSMarine == "1")$OTU
OSM <- subset(otu_overlaps, OSMarine == "1")$OTU
SMO <- subset(otu_overlaps, SMOutfall == "1")$OTU

myV <- plotVenn(list(REstuary1=RE1, REstuary2=RE2, NSMarine=NSM, OSMarine=OSM, SMOutfall=SMO), nCycles = 7000, 
                setColors=c('red', "orange", 'yellow', "blue", "green"), opacity=0.2, labelRegions = F, showNumbers = T, systemShow=TRUE)


#################################################
#Figure 6: Plot relative abundance of shared ASVs
#################################################

###MICROBIAL INDICATORS

#Select only ASVs present at all locations
Mic_pruned.ra.long.shared <- Mic_pruned.ra.long[Mic_pruned.ra.long$OTU %in% rownames(pres_abs_matrix(y.mic))[rowSums(pres_abs_matrix(y.mic))==5],]

#Aggregate by taxonomy
Mic_pruned.ra.long.shared.agg <- aggregate(Abundance~Location+Season+Depth+Class+Family, Mic_pruned.ra.long.shared, FUN= sum)
Mic_pruned.ra.long.shared.agg$Abundance <- Mic_pruned.ra.long.shared.agg$Abundance*100

#Plot
Mic_pruned.ra.long.shared.agg$Location <- factor(Mic_pruned.ra.long.shared.agg$Location, levels = c('R-Estuary-1', 'R-Estuary-2', 'NS-Marine', 'OS-Marine','SM-Outfall'))
Mic_pruned.ra.long.shared.agg$Season<- factor(Mic_pruned.ra.long.shared.agg$Season, levels = c('winter','spring','summer','autumn'))
Mic_pruned.ra.long.shared.agg$Depth<- factor(Mic_pruned.ra.long.shared.agg$Depth, levels = c('surface', 'bottom'))


Mic_shared.otu.plot <- ggplot(Mic_pruned.ra.long.shared.agg, aes(x = Location, y = Abundance, fill = Family)) + 
  facet_grid(Depth~Season, space= "fixed") +
  geom_col()+
  theme(legend.position = "bottom")+ 
  theme_bw()+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportion of shared ASVs") +
  scale_fill_manual(values=indicators.col)+
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1))
Mic_shared.otu.plot

ggsave("./output_graphs/MicIndicatorsSharedASVs.pdf", Mic_shared.otu.plot)


#################################################
#Supplementary: Venn diagram, Mic. indicators
#################################################

#Mic indicators
mi.RE1 <- subset(otu_overlaps.mic, REstuary1 == "1")$OTU
mi.RE2 <- subset(otu_overlaps.mic, REstuary2 == "1")$OTU
mi.NSM <- subset(otu_overlaps.mic, NSMarine == "1")$OTU
mi.OSM <- subset(otu_overlaps.mic, OSMarine == "1")$OTU
mi.SMO <- subset(otu_overlaps.mic, SMOutfall == "1")$OTU

myV <- plotVenn(list(REstuary1=mi.RE1, REstuary2=mi.RE2, NSMarine=mi.NSM, OSMarine=mi.OSM, SMOutfall=mi.SMO), nCycles = 7000, 
                setColors=c('red', "orange", 'yellow', "blue", "green"), opacity=0.2, labelRegions = F, showNumbers = T, systemShow=TRUE)


# CLEAN UP #################################################

# Clear plots
dev.off() 

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
