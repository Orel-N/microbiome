#################################
#Diversity analysis:Dataset 2020
#################################

#Load libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(RColorBrewer)
library(DESeq2)
library(vegan)
library(ggpubr)
library(stats)


#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


############
#Import data
############

#Import data
phy_obj3_20 <- readRDS("./data/phyloseqPrevFiltered20.RDS")
phy_obj3_20


#################
#Plot nMDS 
#################

# DESeq conversion 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
phy_obj3.dds <- phyloseq_to_deseq2(phy_obj3_20, ~1)
geoMeans = apply(counts(phy_obj3.dds), 1, gm_mean)

# DESeq2 Variance Stabilization
phy_obj3.dds = estimateSizeFactors(phy_obj3.dds, geoMeans = geoMeans)
phy_obj3.dds <- estimateDispersions(phy_obj3.dds)
otu.vst <- getVarianceStabilizedData(phy_obj3.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(phy_obj3_20))
ps.vst<-phy_obj3_20
otu_table(ps.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

# Plot the ordination with variance stabilized data
phy_obj3.vst.nmds <- ordinate(ps.vst, method = "NMDS", distance = "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.nmds, label ="Location", color = "Season") + geom_point(aes(shape=Depth), size=3) + 
  stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95) + theme_bw() + scale_color_manual(values=season.col)

ggsave("./output_graphs/nMDS_var.pdf")
saveRDS(ps.vst, "./data/phyloseqVarStab20.RDS")


###########
#PERMANOVA
###########

###statistical significance of the groups
df <- as(sample_data(ps.vst), "data.frame")
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Depth, data= df)
adonis_all

#Post-hoc test 
#posthoc to check which seasons are different
groups <- df[["Season"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

groups <- df[["Location"]]
mod <- betadisper(d, groups)
permutest(mod)

T_Location <- as.data.frame(mod.HSD$group)
write.csv(T_Location, "./output_tables/t_Location.csv")


#############################################
#Preparation of data for taxonomic comparison
#############################################

#Abundance value transformation
phy_obj3ra_20 = transform_sample_counts(phy_obj3_20, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra_20)


########################################################
#Explore community composition on phylum and class level
########################################################

#Calculate abundance for each Phylum
phy_obj3_melt.agg.phylum <- aggregate (Abundance~Location+Season+Depth+Phylum, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.phylum$Abundance <- phy_obj3_melt.agg.phylum$Abundance*100
phy_obj3_melt.agg.phylum<- phy_obj3_melt.agg.phylum[phy_obj3_melt.agg.phylum$Abundance>0,]

phy_obj3_melt.agg.phylum_average <- aggregate (Abundance~Phylum, phy_obj3_melt.agg.phylum, FUN="mean")


#Calculate abundance for each Class
phy_obj3_melt.agg.class <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.class$Abundance <- phy_obj3_melt.agg.class$Abundance*100
phy_obj3_melt.agg.class<- phy_obj3_melt.agg.class[phy_obj3_melt.agg.class$Abundance>0,]

phy_obj3_melt.agg.class.allSamples <- as.data.frame(as.list(aggregate (Abundance~Class, phy_obj3_melt.agg.class, 
                                                 FUN = function(x) c(mean=mean(x), stdev=sd(x)))))
phy_obj3_melt.agg.class.allSeasons <- as.data.frame(as.list(aggregate (Abundance~Season+Class, phy_obj3_melt.agg.class, 
                                                                       FUN = function(x) c(mean=mean(x), stdev=sd(x)))))


#####################################################
#Figure 2: Plot community composition on class level
#####################################################

#Edit order
phy_obj3_melt.agg.class$Season = factor(phy_obj3_melt.agg.class$Season, levels=c( 'winter','spring','summer','autumn'))
phy_obj3_melt.agg.class$Location = factor(phy_obj3_melt.agg.class$Location, levels=c('R-Estuary-1', 'R-Estuary-2', 'NS-Marine', 'OS-Marine','SM-Outfall'))
phy_obj3_melt.agg.class$Depth = factor(phy_obj3_melt.agg.class$Depth, levels=c('surface', 'bottom'))


#remove below 3% ra
threshold<- 3
phy_obj3_melt.agg.class$Class <- as.character(phy_obj3_melt.agg.class$Class)
taxa_classes <- unique(phy_obj3_melt.agg.class$Class[!phy_obj3_melt.agg.class$Abundance<threshold])
phy_obj3_melt.agg.class$Class[phy_obj3_melt.agg.class$Abundance<threshold] <- "Other taxa"
phy_obj3_melt.agg.class$Class <- factor(phy_obj3_melt.agg.class$Class,
                                        levels=c(taxa_classes,"Other taxa"))

#Plot Class(Date)
ggplot(phy_obj3_melt.agg.class, aes(x = Location, y = Abundance, fill = Class))+
  facet_grid(Depth~Season, space= "fixed")+
  theme_bw()+
  scale_fill_manual(values=phyla.col)+
  geom_bar(stat = "identity",position="fill")+
  ylab("Sequence proportion")+
  theme(axis.text.x = element_text(angle = 90))

ggsave("./output_graphs/CommunityCompositionClass.pdf", last_plot())



###############################################################
#Explore community composition on order, family and genus level
###############################################################

#Calculate abundance for each Order
phy_obj3_melt.agg.order <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order, phy_obj3ra.melt,
                                                           FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.order$Abundance <- phy_obj3_melt.agg.order$Abundance.sum*100
phy_obj3_melt.agg.order<- phy_obj3_melt.agg.order[phy_obj3_melt.agg.order$Abundance>0,]
phy_obj3_melt.agg.order.allSamples <- as.data.frame(as.list(aggregate (Abundance~Class+Order, phy_obj3_melt.agg.order, 
                                                                       FUN = function(x) c(mean=mean(x), stdev=sd(x)))))
phy_obj3_melt.agg.order.allSeasons <- as.data.frame(as.list(aggregate (Abundance~Season+Class+Order, phy_obj3_melt.agg.order, 
                                                                       FUN = function(x) c(mean=mean(x), stdev=sd(x)))))


##calculate abundance for each Family
phy_obj3_melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, phy_obj3ra.melt,
                                                            FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.family$Abundance <- phy_obj3_melt.agg.family$Abundance.sum*100
phy_obj3_melt.agg.family<- phy_obj3_melt.agg.family[phy_obj3_melt.agg.family$Abundance>0,]




#Calculate abundance for each Genus
phy_obj3_melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family+Genus, phy_obj3ra.melt,
                                                           FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.genus$Abundance <- phy_obj3_melt.agg.genus$Abundance.sum*100
phy_obj3_melt.agg.genus<- phy_obj3_melt.agg.genus[phy_obj3_melt.agg.genus$Abundance.sum>0,]


#############################################################################
#Supplementary: Plot community composition on family level - selected groups
#############################################################################

#Edit order
phy_obj3_melt.agg.order$Season = factor(phy_obj3_melt.agg.order$Season, levels=c( 'winter','spring','summer','autumn'))
phy_obj3_melt.agg.order$Location = factor(phy_obj3_melt.agg.order$Location, levels=c('R-Estuary-1', 'R-Estuary-2', 'NS-Marine', 'OS-Marine','SM-Outfall'))
phy_obj3_melt.agg.order$Depth = factor(phy_obj3_melt.agg.order$Depth, levels=c('surface', 'bottom'))

###Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj3_melt.agg.order , subset = Class == "Alphaproteobacteria")

#remove below 1% ra
threshold<- 1
Alpha.melt$Order <- as.character(Alpha.melt$Order)
taxa_classes <- unique(Alpha.melt$Order[!Alpha.melt$Abundance<threshold])
Alpha.melt$Order[Alpha.melt$Abundance<threshold] <- "Other taxa"
Alpha.melt$Order <- factor(Alpha.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

colourCount = length(unique(Alpha.melt$Order))
getPalette = colorRampPalette(brewer.pal(11, "Set2"))

ggplot(Alpha.melt, aes(x = Location, y = Abundance, fill = Order))+
  facet_grid(Depth~Season, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  ylab("Sequence proportion")+
  theme(axis.text.x = element_text(angle = 90))

ggsave("./output_graphs/CommunityCompositionAlphaproteobacteria.pdf", last_plot())

###Bacteroidia
Bacter.melt <- subset(phy_obj3_melt.agg.order , subset = Class == "Bacteroidia")

ggplot(Bacter.melt, aes(x = Location, y = Abundance, fill = Order))+
  facet_grid(Depth~Season, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  ylab("Sequence proportion")+
  theme(axis.text.x = element_text(angle = 90))  

ggsave("./output_graphs/CommunityCompositionBacteroidia.pdf", last_plot())

### Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj3_melt.agg.order , subset = Class == "Gammaproteobacteria")

#remove below 1% ra
threshold<- 1
Gamma.melt$Order <- as.character(Gamma.melt$Order)
taxa_classes <- unique(Gamma.melt$Order[!Gamma.melt$Abundance<threshold])
Gamma.melt$Order[Gamma.melt$Abundance<threshold] <- "Other taxa"
Gamma.melt$Order <- factor(Gamma.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

colourCount = length(unique(Gamma.melt$Order))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

ggplot(Gamma.melt, aes(x = Location, y = Abundance, fill = Order))+
  facet_grid(Depth~Season, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  ylab("Sequence proportion")+
  theme(axis.text.x = element_text(angle = 90)) 

ggsave("./output_graphs/CommunityCompositionGammaproteobacteria.pdf", last_plot())



###############################
#Microbial pollution indicators
###############################

#Select microbial indicators
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra_20, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)


#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), asv_number=length(which(x!=0))))))
#Explore dataset
length(unique(ps_Mic_Ind_ra.melt$OTU))
ps_Mic_Ind_ra.melt.agg.sampling <- as.data.frame(as.list(aggregate (cbind(Abundance.sum, Abundance.asv_number)~Location+Season+Depth, ps_Mic_Ind_ra.melt.agg.family, 
                                                                    FUN=function(x) c(sum=sum(x), families=length(which(x!=0))))))
ps_Mic_Ind_ra.melt.agg.location.average <- as.data.frame(as.list(aggregate (cbind(Abundance.sum, Abundance.asv_number)~Location, ps_Mic_Ind_ra.melt.agg.sampling, 
                                                                            FUN=function(x) c(mean = mean(x), stdev=sd(x)))))
ps_Mic_Ind_ra.melt.agg.season.average <- as.data.frame(as.list(aggregate (cbind(Abundance.sum, Abundance.asv_number)~Season, ps_Mic_Ind_ra.melt.agg.sampling, 
                                                                          FUN=function(x) c(mean = mean(x), stdev=sd(x)))))

######################################
# Figure 3: Heatmap of mic. indicators
######################################

library(rstatix)

ps_Mic_Ind_ra.melt.agg.family$Abundance.perc <- ps_Mic_Ind_ra.melt.agg.family$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.family$Abundance.sum <- format(ps_Mic_Ind_ra.melt.agg.family$Abundance.sum, scientific = FALSE)

ps_Mic_Ind_ra.melt.agg.family$Location = factor(ps_Mic_Ind_ra.melt.agg.family$Location, levels=c('R-Estuary-1', 'R-Estuary-2', 'NS-Marine', 'OS-Marine','SM-Outfall'))
ps_Mic_Ind_ra.melt.agg.family$Depth = factor(ps_Mic_Ind_ra.melt.agg.family$Depth, levels=c('surface', 'bottom'))
ps_Mic_Ind_ra.melt.agg.family$Season = factor(ps_Mic_Ind_ra.melt.agg.family$Season, levels=c('winter', 'spring', 'summer', 'autumn'))
ps_Mic_Ind_ra.melt.agg.family$Abundance.class <- ifelse(ps_Mic_Ind_ra.melt.agg.family$Abundance.perc == "0", "0", ifelse(ps_Mic_Ind_ra.melt.agg.family$Abundance.sum < "0.005", "1", 
                                                        ifelse(ps_Mic_Ind_ra.melt.agg.family$Abundance.sum < "0.01", "2", 
                                                               ifelse(ps_Mic_Ind_ra.melt.agg.family$Abundance.sum < "0.05", "3", "4"))))
                                                                            # ifelse(ps_Mic_Ind_ra.melt.agg.family$Abundance.perc < "1", "5","6"))))))
ps_Mic_Ind_ra.melt.agg.family$Abundance.class <- as.numeric(ps_Mic_Ind_ra.melt.agg.family$Abundance.class)

heatmap <- ggplot(ps_Mic_Ind_ra.melt.agg.family, aes(Season, Family, fill= Abundance.class)) + 
  geom_tile() +
  xlab(label = "Sample") +
  theme_bw()+
  facet_grid (Depth~Location, switch = "both", scales = "free", space = "free") +
  scale_fill_gradient2(name = "Sequence proportion",
                       high = "red",
                       mid = "white",
                       na.value = "white",
                       labels = c("0", "<0.5%", "0.5%-1%", "1%-5%", ">5%")) +
  theme(strip.text.y = element_text(angle=0),
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        panel.grid = element_blank())
heatmap

ggsave("./output_graphs/MicrobialIndicatorsHeatmap.pdf", last_plot())


#################################
#Mic. indicators in both datasets
#################################

#Import data
phy_obj3 <- readRDS("./data/phyloseqPrevFiltered.RDS")
phy_obj3

#Abundance value transformation
phy_obj3ra = transform_sample_counts(phy_obj3, function (otu) {otu/sum(otu)})

ps_Mic_Ind_15_20_ra <- subset_taxa(phy_obj3ra, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_15_20_ra.melt <- psmelt(ps_Mic_Ind_15_20_ra)


#Calculate abundance for each taxa
ps_Mic_Ind_15_20_ra.melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Dataset+Class+Order+Family, ps_Mic_Ind_15_20_ra.melt,
                                                                 FUN = function(x) c(sum = sum(x), asv_number=length(which(x!=0))))))

#Plot
ps_Mic_Ind_15_20_ra.melt$Season = factor(ps_Mic_Ind_15_20_ra.melt$Season, levels=c('winter','spring','summer','autumn'))
ps_Mic_Ind_15_20_ra.melt$Location = factor(ps_Mic_Ind_15_20_ra.melt$Location, levels=c('R-Mouth', 'R-Estuary-1', 'R-Estuary-2', 'NS-Marine', 'OS-Marine','SM-Outfall'))
ps_Mic_Ind_15_20_ra.melt$Date = factor(ps_Mic_Ind_15_20_ra.melt$Date, levels=c('4.12.2018','25.04.2019','18.06.2019','15.10.2019','4.02.2016', '24.03.2015', '15.07.2015', '25.11.2015'))
ps_Mic_Ind_15_20_ra.melt$Depth = factor(ps_Mic_Ind_15_20_ra.melt$Depth, levels = c('surface', 'bottom'))

ggplot(ps_Mic_Ind_15_20_ra.melt, aes(x = Location, y = Abundance, fill = Family))+
  facet_grid(Depth~Date, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=indicators.col)+
  ylab("Sequence proportion")+
  theme(axis.text.x = element_text(angle = 90))

ggsave("./output_graphs/MicrobialIndicators_15_20.pdf", last_plot())

#Explore 2015 dataset
ps_Mic_Ind_15_ra.melt <- filter(ps_Mic_Ind_15_20_ra.melt, Dataset == "2015")
length(unique(ps_Mic_Ind_15_ra.melt$OTU))

ps_Mic_Ind_15_ra.melt.agg.family <- filter(ps_Mic_Ind_15_20_ra.melt.agg.family, Dataset == "2015")
ps_Mic_Ind_15_ra.melt.agg.sampling <- as.data.frame(as.list(aggregate (cbind(Abundance.sum, Abundance.asv_number)~Location+Season+Depth, ps_Mic_Ind_15_ra.melt.agg.family, 
                                                                    FUN=function(x) c(sum=sum(x), families=length(which(x!=0))))))

#####################
#Arcobacter phylogeny
#####################
library(Biostrings)


###Subset taxa and extract seqences

##Only 2020
Arcobacter_family_20 <- subset_taxa(phy_obj3ra_20, Family == "Arcobacteraceae")
saveRDS(Arcobacter_family_20, "./data/Arcobacter_family_ps_20.RDS")

Arcobacter_family.rs_20 <- refseq(Arcobacter_family_20)
writeXStringSet(Arcobacter_family.rs_20, "./data/Arcobacter_family_20.fasta", format = "fasta")

##Only 2015
ps_Mic_Ind_15_ra <- subset_samples(phy_obj3ra, Dataset == "2015")

Arcobacter_family_15 <- subset_taxa(ps_Mic_Ind_15_ra, Family == "Arcobacteraceae")
saveRDS(Arcobacter_family_15, "./data/Arcobacter_family_ps_15.RDS")

##Both datsets
Arcobacter_family_15_20 <- subset_taxa(ps_Mic_Ind_15_20_ra, Family == "Arcobacteraceae")
saveRDS(Arcobacter_family_15_20, "./data/Arcobacter_family_ps_15_20.RDS")

Arcobacter_family.rs_15_20 <- refseq(Arcobacter_family_15_20)
writeXStringSet(Arcobacter_family.rs_15_20, "./data/Arcobacter_family_15_20.fasta", format = "fasta")

#Explore dataset
Arcobacter <- filter(ps_Mic_Ind_15_20_ra.melt, Family == "Arcobacteraceae")


ASV55 <- filter(ps_Mic_Ind_15_20_ra.melt, OTU == "ASV55")


ggplot(data=subset(ps_Mic_Ind_15_20_ra.melt, OTU == "ASV55"), aes(x = Abundance*100, y = Location, fill = as.factor(Dataset)))+
  facet_grid(Date~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  xlab("Relative seq. abundance [%]")+
  labs(fill = "Dataset")


###################################################
#Figure 6: Phylogenetic tree with closest relatives
###################################################

#Upload Arcobacter_family.fasta on Silva webside and create tree 
#(Classify 0.97 similarity, compute tree: Denovo including neighbours, program:RAxML, model: GTR)
#View in wasabi and download results
#Export tree in Newick format
#Visualize using iTOL



# CLEAN UP #################################################

# Clear plots
dev.off() 

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
