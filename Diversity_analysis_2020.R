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

#Edit data - subset based on dataset
phy_obj3_20 <- subset_samples(phy_obj3, Dataset == "2020")

####################################
#Alpha diversity plot
####################################

#Create Alphadiversity table 2020
alpha_df <- estimate_richness(phy_obj3_20, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
alpha_df_20 <- data.frame(sample_data(phy_obj3_20), alpha_df)
write.table(alpha_df_20, "./tables/AlphaDiversity_20.txt")


#Edit table 2020
alpha_df_20$Season <- factor(alpha_df_20$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_20$Location <- factor(alpha_df_20$Location, levels = c("R-Mouth", "R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_20$Depth <- factor(alpha_df_20$Depth, levels = c("bottom", "surface"))
alpha_df_20$Pollution <- ifelse(alpha_df_20$Location == "OS-Marine", "Unpolluted", 
                                   ifelse(alpha_df_20$Location == "NS-Marine", "Unpolluted", "Polluted"))

#Plot Chao1 - separate for depth 2020
Chao.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.20.p
ggsave("./output_graphs/Chao_20.pdf")

#Plot Shannon d.i. - separate for depth 2020
Shannon.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.20.p
ggsave("./output_graphs/Shannon_20.pdf")

#Plot InvSimpson - separate for depth 2020
InvSimpson.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=InvSimpson, color=Location)) +
  geom_point(mapping=aes(x=Location, y=InvSimpson, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
InvSimpson.20.p
ggsave("./output_graphs/InvSimpsin_20.pdf")

alpha.20 <- ggarrange(Chao.20.p, Shannon.20.p, InvSimpson.20.p, nrow=1, common.legend = TRUE)
alpha.20

###Statistic 2020
#Subsetting data 2020 - bottom vs.surface
alpha_df_20b <- filter(alpha_df_20, Depth == "bottom")
alpha_df_20s <- filter(alpha_df_20, Depth == "surface")

#Shapiro-Wilk's test
hist(alpha_df_20$Chao1)
qqnorm(alpha_df_20$Chao1)
qqline(alpha_df_20$Chao1)
hist(alpha_df_20$Shannon)
qqnorm(alpha_df_20$Shannon)
qqline(alpha_df_20$Shannon)

alpha_df_20 %>% shapiro_test(Chao1, Shannon, InvSimpson)
alpha_df_20b %>% shapiro_test(Chao1, Shannon, InvSimpson)
alpha_df_20s %>% shapiro_test(Chao1, Shannon, InvSimpson)

#ANOVA or Kruskal-Wallis test
anova.Ch.20 <- aov(alpha_df_20$Chao1 ~ Depth+Season+Pollution+Location, data=alpha_df_20)
summary(anova.Ch.20)
TukeyHSD(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20b$Chao1 ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20s$Chao1 ~ Season+Pollution+Location, data=alpha_df_20s)
summary(anova.Ch.20)

kruskal.test(alpha_df_20$Shannon ~ Season, data = alpha_df_20) 
kruskal.test(alpha_df_20$Shannon ~ Depth, data = alpha_df_20) 
kruskal.test(alpha_df_20$Shannon ~ Location, data = alpha_df_20)
kruskal.test(alpha_df_20$Shannon ~ Pollution, data = alpha_df_20)
anova.Ch.20 <- aov(alpha_df_20b$Shannon ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20s$Shannon ~ Season+Pollution+Location, data=alpha_df_20s)
summary(anova.Ch.20)
TukeyHSD(anova.Ch.20)

kruskal.test(alpha_df_20$InvSimpson ~ Season, data = alpha_df_20) 
kruskal.test(alpha_df_20$InvSimpson ~ Location, data = alpha_df_20)
kruskal.test(alpha_df_20$InvSimpson ~ Pollution, data = alpha_df_20)
anova.InvS.20 <- aov(alpha_df_20$InvSimpson ~ Season+Pollution+Location, data=alpha_df_20)
summary(anova.InvS.20)

anova.InvS.20b <- aov(alpha_df_20b$InvSimpson ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.InvS.20b)

kruskal.test(alpha_df_20s$InvSimpson ~ Season, data = alpha_df_20s) 
kruskal.test(alpha_df_20s$InvSimpson ~ Location, data = alpha_df_20s)
kruskal.test(alpha_df_20s$InvSimpson ~ Pollution, data = alpha_df_20s)


##############################
#Preparation of data for taxonomic comparison
##############################

#Abundance value transformation
phy_obj3ra_20 = transform_sample_counts(phy_obj3_20, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra_20.melt <- psmelt(phy_obj3ra_20)
write.table(phy_obj3ra_20.melt, "./tables/phy_obj3ra.melt.20.txt")


###########################
#Taxonomic compositions on class level 2020 only!
###########################

phy_obj3ra.melt <- phy_obj3ra_20.melt

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
colourCount = length(unique(phy_obj3_melt.agg.class$Class))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
phy_obj3_melt.agg.class$Season = factor(phy_obj3_melt.agg.class$Season, levels=c( 'winter','spring','summer','autumn'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))


###Subsets by taxonomy - Campylobacteria (position stack/position fill)
Campyl.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Campylobacteria")

Campyl.melt$Season_f = factor(Campyl.melt$Season, levels=c( 'winter','spring','summer','autumn'))
ggplot(Campyl.melt, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Season_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
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
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
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
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
ggsave("./output_graphs/Gammaproteobacteria.pdf", last_plot())


###Subsets by taxonomy - Bacteroidia
Bacter.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacteroidia")
Bacter.melt$Season_f = factor(Bacter.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Bacter.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Bacilli
Bacilli.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacilli")
Bacilli.melt$Season_f = factor(Bacilli.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Bacilli.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Clostridia
Clost.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Clostridia")
Clost.melt$Season_f = factor(Clost.melt$Season, levels=c('winter','spring','summer','autumn'))

ggplot(Clost.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
ggsave("./output_graphs/Clostridia.pdf", last_plot())

###Vibrio
Vibrio<-subset(phy_obj3_melt.agg.genus, subset = Genus == "Vibrio")

##################################
#Microbial pollution indicators
##################################

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra_20, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)
write.table(ps_Mic_Ind_ra.melt, "./tables/ps_Mic_Ind_ra.melt.txt")

#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Sample+Season+Date+Class+Order+Family+Genus, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]

ps_Mic_Ind_ra.melt.agg.genus$Season_f = factor(ps_Mic_Ind_ra.melt.agg.genus$Season, levels=c('winter','spring','summer','autumn'))
ps_Mic_Ind_ra.melt.agg.genus$Location_f = factor(ps_Mic_Ind_ra.melt.agg.genus$Location, levels=c('R-Estuary-1','R-Estuary-2', 'NS-Marine', 'OS-Marine', 'SM-Outfall'))

colourCount = length(unique(ps_Mic_Ind_ra.melt.agg.genus$Family))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
ggplot(ps_Mic_Ind_ra.melt.agg.genus, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Season_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1'))
ggsave("./output_graphs/MicrobialInidicators.pdf", last_plot())

#Separated graphs
Mic_Ind_1 <- subset(ps_Mic_Ind_ra.melt.agg.genus, (Family %in% c("Arcobacteraceae", "Aeromonadaceae", "Pseudomonadaceae")))
Mic_Ind_2 <- subset(ps_Mic_Ind_ra.melt.agg.genus, (Family %in% c("Clostridiaceae", "Enterobacteriaceae")))
Mic_Ind_3 <- subset(ps_Mic_Ind_ra.melt.agg.genus, (Family %in% c("Legionellaceae", "Leptospiraceae", "Mycobacteriaceae")))
Mic_Ind_4 <- subset(ps_Mic_Ind_ra.melt.agg.genus, (Family %in% c("Chlamydiaceae", "Moraxellaceae")))
Mic_Ind_5 <- subset(ps_Mic_Ind_ra.melt.agg.genus, (Family %in% c("Enterobacteriaceae", "Enterococcaceae")))


colourCount = length(unique(ps_Mic_Ind_ra.melt.agg.genus$Family))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
ggplot(Mic_Ind_1, aes(x = Season_f, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Depth~Family+Location_f, strip.position = "bottom", nrow = 2, drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
  

s1 <- ggplot(Mic_Ind_1[Mic_Ind_1$Depth=="surface",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
s1

b1 <- ggplot(Mic_Ind_1[Mic_Ind_1$Depth=="bottom",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
b1

ggarrange(s1, b1, nrow=2)

s2 <- ggplot(Mic_Ind_2[Mic_Ind_2$Depth=="surface",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
s2

b2 <- ggplot(Mic_Ind_2[Mic_Ind_2$Depth=="bottom",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
b2

ggarrange(s2, b2, nrow=2)

s3 <- ggplot(Mic_Ind_3[Mic_Ind_3$Depth=="surface",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
s3

b3 <- ggplot(Mic_Ind_3[Mic_Ind_3$Depth=="bottom",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
b3

ggarrange(s3, b3, nrow=2)

s4 <- ggplot(Mic_Ind_4[Mic_Ind_4$Depth=="surface",], aes(x = Season_f, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
s4

b4 <- ggplot(Mic_Ind_4[Mic_Ind_4$Depth=="bottom",], aes(x = Season_f, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
b4

ggarrange(s4, b4, nrow=2)

s5 <- ggplot(Mic_Ind_5[Mic_Ind_5$Depth=="surface",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
s5

b5 <- ggplot(Mic_Ind_5[Mic_Ind_5$Depth=="bottom",], aes(x = Season_f, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="stack") +
  facet_wrap(Family~Location_f, nrow=1, strip.position = "bottom", drop=FALSE) +
  theme_classic() +
  theme(strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
b5

ggarrange(s5, b5, nrow=2)

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
adonis_all <- adonis2(d ~ Pollution+Location+Season+Depth, data= df, perm = 999)
adonis_all

#Post-hoc test (??)


df <- as(sample_data(ps.ind.vst), "data.frame")
df$Pollution <- ifelse(df$Location == "OS-Marine", "Unpolluted", 
                       ifelse(df$Location == "NS-Marine", "Unpolluted", "Polluted"))
d <- phyloseq::distance(ps.ind.vst, "euclidean")
adonis_all <- adonis2(d ~ Pollution+Location+Season+Depth+Dataset, data= df, perm = 999)
adonis_all
#Post-hoc test (??)

