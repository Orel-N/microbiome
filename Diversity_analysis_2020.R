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

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


############
#Import data
############

#Import data
phy_obj3 <- readRDS("./data/phyloseqPrevFiltered.RDS")
phy_obj3

#Edit data - subset based on dataset
phy_obj3_20 <- subset_samples(phy_obj3, Dataset == "2020")
phy_obj3_20 <- prune_taxa(taxa_sums(phy_obj3_20)>0, phy_obj3_20)
phy_obj3_20


#############################################
#Preparation of data for taxonomic comparison
#############################################

#Abundance value transformation
phy_obj3ra_20 = transform_sample_counts(phy_obj3_20, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra_20)


################################################################
#Taxonomic compositions on phylum, class, family and genus level
################################################################

#Calculate abundance for each Phylum
phy_obj3_melt.agg.phylum <- aggregate (Abundance~Location+Season+Depth+Phylum, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.phylum$Abundance <- phy_obj3_melt.agg.phylum$Abundance*100
phy_obj3_melt.agg.phylum<- phy_obj3_melt.agg.phylum[phy_obj3_melt.agg.phylum$Abundance>0,]

phy_obj3_melt.agg.phylum_average <- aggregate (Abundance~Phylum, phy_obj3_melt.agg.phylum, FUN="mean")

##calculate abundance for each Family
phy_obj3_melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, phy_obj3ra.melt,
                                                            FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.family$Abundance <- phy_obj3_melt.agg.family$Abundance.sum*100
phy_obj3_melt.agg.family<- phy_obj3_melt.agg.family[phy_obj3_melt.agg.family$Abundance>0,]

#Calculate abundance for each taxa
phy_obj3_melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family+Genus, phy_obj3ra.melt,
                                                           FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.genus$Abundance <- phy_obj3_melt.agg.genus$Abundance.sum*100
phy_obj3_melt.agg.genus<- phy_obj3_melt.agg.genus[phy_obj3_melt.agg.genus$Abundance.sum>0,]

#Calculate abundance for each Class
phy_obj3_melt.agg.class <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.class$Abundance <- phy_obj3_melt.agg.class$Abundance*100
phy_obj3_melt.agg.class<- phy_obj3_melt.agg.class[phy_obj3_melt.agg.class$Abundance>0,]

phy_obj3_melt.agg.class.allSamples <- aggregate (Abundance~Class, phy_obj3_melt.agg.class, FUN="mean")


#Plot community composition on class level
##########################################

#Edit order
phy_obj3_melt.agg.class$Season = factor(phy_obj3_melt.agg.class$Season, levels=c( 'winter','spring','summer','autumn'))
phy_obj3_melt.agg.class$Location = factor(phy_obj3_melt.agg.class$Location, levels=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))

#remove below 3% ra
threshold<- 3
phy_obj3_melt.agg.class$Class <- as.character(phy_obj3_melt.agg.class$Class)
taxa_classes <- unique(phy_obj3_melt.agg.class$Class[!phy_obj3_melt.agg.class$Abundance<threshold])
phy_obj3_melt.agg.class$Class[phy_obj3_melt.agg.class$Abundance<threshold] <- "Other taxa"
phy_obj3_melt.agg.class$Class <- factor(phy_obj3_melt.agg.class$Class,
                                        levels=c(taxa_classes,"Other taxa"))

#Plot Class(Date)
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  theme_bw()+
  scale_fill_manual(values=phyla.col)+
  geom_bar(stat = "identity",position="fill")

ggsave("./output_graphs/CommunityCompositionClass.pdf", last_plot())


#Plot community composition on family level - selected groups
#############################################################

#Edit order
phy_obj3_melt.agg.family$Season = factor(phy_obj3_melt.agg.family$Season, levels=c( 'winter','spring','summer','autumn'))
phy_obj3_melt.agg.family$Location = factor(phy_obj3_melt.agg.family$Location, levels=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))

###Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Alphaproteobacteria")

#remove below 1% ra
threshold<- 1
Alpha.melt$Order <- as.character(Alpha.melt$Order)
taxa_classes <- unique(Alpha.melt$Order[!Alpha.melt$Abundance<threshold])
Alpha.melt$Order[Alpha.melt$Abundance<threshold] <- "Other taxa"
Alpha.melt$Order <- factor(Alpha.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

colourCount = length(unique(Alpha.melt$Order))
getPalette = colorRampPalette(brewer.pal(11, "Set2"))

ggplot(Alpha.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))

ggsave("./output_graphs/CommunityCompositionAlphaproteobacteria.pdf", last_plot())

### Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Gammaproteobacteria")

#remove below 1% ra
threshold<- 1
Gamma.melt$Order <- as.character(Gamma.melt$Order)
taxa_classes <- unique(Gamma.melt$Order[!Gamma.melt$Abundance<threshold])
Gamma.melt$Order[Gamma.melt$Abundance<threshold] <- "Other taxa"
Gamma.melt$Order <- factor(Gamma.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

colourCount = length(unique(Gamma.melt$Order))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

ggplot(Gamma.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))

ggsave("./output_graphs/CommunityCompositionGammaproteobacteria.pdf", last_plot())

###Bacteroidia
Bacter.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacteroidia")

ggplot(Bacter.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")

ggsave("./output_graphs/CommunityCompositionBacteroidia.pdf", last_plot())


###############################
#Microbial pollution indicators
###############################

#Select microbial indicators
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra_20, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)

#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]

ps_Mic_Ind_ra.melt.agg.genus$Season = factor(ps_Mic_Ind_ra.melt.agg.genus$Season, levels=c('winter','spring','summer','autumn'))
ps_Mic_Ind_ra.melt.agg.genus$Location = factor(ps_Mic_Ind_ra.melt.agg.genus$Location, levels=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))

ggplot(ps_Mic_Ind_ra.melt.agg.genus, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=indicators.col)

ggsave("./output_graphs/MicrobialIndicators.pdf", last_plot())

# Heatmap
ps_Mic_Ind_ra_glom = tax_glom(ps_Mic_Ind_ra, taxrank="Family")
fig <- plot_heatmap(ps_Mic_Ind_ra_glom, "Family", taxa.label = "Family", sample.label="Location", taxa.order="Class", sample.order="Location")
fig


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
saveRDS(ps.vst, "./data/phyloseqVarStab_20.RDS")

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


# CLEAN UP #################################################

# Clear plots
dev.off() 

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
