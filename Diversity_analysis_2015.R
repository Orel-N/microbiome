#MICROBIOME - Diversity analysis separate for 2015 and 2020 dataset
#Neža Orel, neza.orel@nib.si
#Script for alpha and beta diversity analyses: graphs and statistic
#Data: project MICROBIOME, dataset 2015 and 2020 seperatelly

############################################################################

#Set working directory
setwd("C:/Users/nezao/Documents/5-R/Microbiome2015_2020")

#Load packages
library(ggplot2)
library(dplyr)
library(ggpubr)
library(phyloseq)
library(RColorBrewer)
library(DESeq2)
library(rstatix)
library(vegan)


#Import data
phy_obj3 <- readRDS("./phyloseqFiltered.RDS")
phy_obj3

#prevdf3 <- readRDS("./prevdfFiltered.RDS")
#summary(prevdf3)
#str(prevdf3)

#Edit data - subset based on dataset
phy_obj3_15 <- subset_samples(phy_obj3, Dataset == "2015")

#Define colors
phyla.col <- readRDS("./data/phyla_col.RDS")

indicators.col <- readRDS("./data/indicators_col.RDS")

####################################
#Alpha diversity plot
#

##############################
#Preparation of data for taxonomic comparison
##############################

#Abundance value transformation
phy_obj3ra_15 = transform_sample_counts(phy_obj3_15, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra_15.melt <- psmelt(phy_obj3ra_15)
write.table(phy_obj3ra_15.melt, "./tables/phy_obj3ra.melt.15.txt")



###########################
#Taxonomic compositions on class level 2020 only!
###########################

phy_obj3ra.melt <- phy_obj3ra_15.melt

#Calculate abundance for each taxa
phy_obj3_melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family+Genus, phy_obj3ra.melt,
                                                           FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.genus$Abundance <- phy_obj3_melt.agg.genus$Abundance.sum*100
phy_obj3_melt.agg.genus<- phy_obj3_melt.agg.genus[phy_obj3_melt.agg.genus$Abundance.sum>0,]


##calculate abundance for each Family
phy_obj3_melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, phy_obj3ra.melt,
                                                            FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_melt.agg.family$Abundance <- phy_obj3_melt.agg.family$Abundance.sum*100
phy_obj3_melt.agg.family<- phy_obj3_melt.agg.family[phy_obj3_melt.agg.family$Abundance.sum>0,]


#Calculate abundance for each Class
phy_obj3_melt.agg.class <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.class$Abundance <- phy_obj3_melt.agg.class$Abundance*100
phy_obj3_melt.agg.class<- phy_obj3_melt.agg.class[phy_obj3_melt.agg.class$Abundance>0,]

#remove below 3% ra
threshold<- 3
phy_obj3_melt.agg.class$Class <- as.character(phy_obj3_melt.agg.class$Class)
taxa_classes <- unique(phy_obj3_melt.agg.class$Class[!phy_obj3_melt.agg.class$Abundance<threshold])
phy_obj3_melt.agg.class$Class[phy_obj3_melt.agg.class$Abundance<threshold] <- "Other taxa"
phy_obj3_melt.agg.class$Class <- factor(phy_obj3_melt.agg.class$Class,
                                        levels=c(taxa_classes,"Other taxa"))
phy_obj3_melt.agg.class.table <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj3_melt.agg.class, FUN="sum")
write.csv(phy_obj3_melt.agg.class.table, "./tables/Class.csv")


#Plot Class(Date)
phy_obj3_melt.agg.class$Season = factor(phy_obj3_melt.agg.class$Season, levels=c( 'winter','spring','summer','autumn'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=phyla.col)+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))


###Subsets by taxonomy - Campylobacteria (position stack/position fill)
Campyl.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Campylobacteria")

Campyl.melt$Season_f = factor(Campyl.melt$Season, levels=c( 'winter','spring','summer','autumn'))
ggplot(Campyl.melt, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Season_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Campylobacteria.pdf", last_plot())


###Subsets by taxonomy - Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Alphaproteobacteria")
Alpha.melt$Season_f = factor(Alpha.melt$Season, levels=c( 'winter','spring','summer','autumn'))

colourCount = length(unique(Alpha.melt$Order))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

#remove below 1% ra
threshold<- 1
Alpha.melt$Order <- as.character(Alpha.melt$Order)
taxa_classes <- unique(Alpha.melt$Order[!Alpha.melt$Abundance<threshold])
Alpha.melt$Order[Alpha.melt$Abundance<threshold] <- "Other taxa"
Alpha.melt$Order <- factor(Alpha.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

ggplot(Alpha.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Alphaproteobacteria.pdf", last_plot())

###Subsets by taxonomy - Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Gammaproteobacteria")
Gamma.melt$Season_f = factor(Gamma.melt$Season, levels=c('winter','spring','summer','autumn'))


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
  facet_grid(Season_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Gammaproteobacteria.pdf", last_plot())


###Subsets by taxonomy - Bacteroidia
Bacter.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacteroidia")
Bacter.melt$Season_f = factor(Bacter.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Bacter.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Bacilli
Bacilli.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacilli")
Bacilli.melt$Season_f = factor(Bacilli.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Bacilli.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Clostridia
Clost.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Clostridia")
Clost.melt$Season_f = factor(Clost.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Clost.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/Clostridia.pdf", last_plot())

###Vibrio
Vibrio<-subset(phy_obj3_melt.agg.genus, subset = Genus == "Vibrio")

##################################
#Microbial pollution indicators
##################################

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra_15, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)
write.table(ps_Mic_Ind_ra.melt, "./tables/ps_Mic_Ind_ra.melt.txt")

#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]

ps_Mic_Ind_ra.melt.agg.genus$Season_f = factor(ps_Mic_Ind_ra.melt.agg.genus$Season, levels=c('winter','spring','summer','autumn'))
colourCount = length(unique(ps_Mic_Ind_ra.melt.agg.genus$Family))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
ggplot(ps_Mic_Ind_ra.melt.agg.genus, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('NS-Marine','R-Estuary-1', 'R-Mouth'))
ggsave("./output_graphs/MicrobialInidicators.pdf", last_plot())

# Heatmap
ps_Mic_Ind_ra_glom = tax_glom(ps_Mic_Ind_ra, taxrank="Family")
fig <- plot_heatmap(ps_Mic_Ind_ra_glom, "Family", taxa.label = "Family", sample.label="Season", taxa.order="Class", sample.order="Location")
fig


#################################
#Plot nMDS / RDA
#################################
# DESeq conversion 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
phy_obj3.dds <- phyloseq_to_deseq2(phy_obj3_15, ~1)
geoMeans = apply(counts(phy_obj3.dds), 1, gm_mean)

# DESeq2 Variance Stabilization
phy_obj3.dds = estimateSizeFactors(phy_obj3.dds, geoMeans = geoMeans)
phy_obj3.dds <- estimateDispersions(phy_obj3.dds)
otu.vst <- getVarianceStabilizedData(phy_obj3.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(phy_obj3_15))
ps.vst<-phy_obj3_15
otu_table(ps.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

# Do the ordination with variance stabilized data
phy_obj3.vst.nmds <- ordinate(ps.vst, "NMDS", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.nmds, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/nMDS_var.pdf")

phy_obj3.vst.rda <- ordinate(ps.vst, "RDA", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.rda, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/RDA_var.pdf")
plot_ordination(ps.vst, phy_obj3.vst.rda, type="taxa", color="Phylum", title="taxa")

#########################
#Microbial indicators

ps.ind.vst <- subset_taxa(ps.vst, Family %in% c(Mic_Ind$Family))

# Do the ordination with variance stabilized data
phy_ind.vst.nmds <- ordinate(ps.ind.vst, "NMDS", "euclidean")
plot_ordination(ps.ind.vst, phy_ind.vst.nmds, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/nMDS_Ind_var.pdf")

phy_ind.vst.rda <- ordinate(ps.ind.vst, "RDA", "euclidean")
plot_ordination(ps.ind.vst, phy_ind.vst.rda, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/RDA_Ind_var.pdf")
plot_ordination(ps.ind.vst, phy_ind.vst.rda, type="taxa", color="Family", title="taxa")+ facet_wrap(~Class, 3)

plot_ordination(ps.ind.vst, phy_ind.vst.rda, type="split", color="Family", title="taxa", label="Location")+ facet_wrap(~Phylum, 3)

###########################
#PERMANOVA

#statistical significance of the groups
df <- as(sample_data(ps.vst), "data.frame")
df$Pollution <- ifelse(df$Location == "OS-Marine", "Unpolluted", 
                       ifelse(df$Location == "NS-Marine", "Unpolluted", "Polluted"))
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Pollution+Location+Season, data= df, perm = 999)
adonis_all

#Post-hoc test (??)


df <- as(sample_data(ps.ind.vst), "data.frame")
df$Pollution <- ifelse(df$Location == "OS-Marine", "Unpolluted", 
                       ifelse(df$Location == "NS-Marine", "Unpolluted", "Polluted"))
d <- phyloseq::distance(ps.ind.vst, "euclidean")
adonis_all <- adonis2(d ~ Pollution+Location+Season, data= df, perm = 999)
adonis_all
#Post-hoc test (??)

