#MICROBIOME - Diversity analysis 2015 and 2020 dataset
#Neža Orel, neza.orel@nib.si
#Script for alpha and beta diversity analyses: graphs and statistic
#Data: project MICROBIOME, dataset 2015 and 2020 together

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

prevdf3 <- readRDS("./prevdfFiltered.RDS")
summary(prevdf3)
str(prevdf3)


####################################
#Alpha diversity plot
####################################

#Create Alphadiversity table
alpha_df <- estimate_richness(phy_obj3, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
alpha_df_15.20 <- data.frame(sample_data(phy_obj3), alpha_df)
write.table(alpha_df_15.20, "./tables/AlphaDiversity15_20.txt")

#Edit table
alpha_df_15.20$Season <- factor(alpha_df_15.20$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_15.20$Location <- factor(alpha_df_15.20$Location, levels = c("R-Mouth", "R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_15.20$Depth <- factor(alpha_df_15.20$Depth, levels = c("bottom", "surface"))
alpha_df_15.20$Pollution <- ifelse(alpha_df_15.20$Location == "OS-Marine", "Unpolluted", 
                             ifelse(alpha_df_15.20$Location == "NS-Marine", "Unpolluted", "Polluted"))
#Plot Chao1 - separate for depth
Chao.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.p
ggsave("./output_graphs/Chao_15_20.pdf")

#Plot Shannon d.i. - separate for depth
Shannon.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.p
ggsave("./output_graphs/Shannon_15_20.pdf")

#Plot InvSimpson - separate for depth
InvSimpson.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=InvSimpson, color=Location)) +
  geom_point(mapping=aes(x=Location, y=InvSimpson, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
InvSimpson.p
ggsave("./output_graphs/InvSimpsin_15_20.pdf")

ggarrange(Chao.p, Shannon.p, InvSimpson.p, nrow=1, common.legend = TRUE)

###Statistic

#Shapiro-Wilk's test
hist(alpha_df_15.20$Chao1)
qqnorm(alpha_df_15.20$Chao1)
qqline(alpha_df_15.20$Chao1)
hist(alpha_df_15.20$Shannon)
qqnorm(alpha_df_15.20$Shannon)
qqline(alpha_df_15.20$Shannon)

alpha_df_15.20 %>% shapiro_test(Chao1, Shannon, InvSimpson)

anova.Ch.15.20 <- aov(alpha_df_15.20$Chao1 ~ Depth+Season+Pollution+Location, data=alpha_df_15.20)
summary(anova.Ch.15.20)

anova.Ch.15.20 <- aov(alpha_df_15.20$Chao1 ~ Location+Depth, data=alpha_df_15.20)
summary(anova.Ch.15.20)
anova.Ch.15.20 <- aov(alpha_df_15.20$Chao1 ~ Location*Depth, data=alpha_df_15.20)
summary(anova.Ch.15.20)
TukeyHSD(anova.Ch.15.20)

kruskal.test(alpha_df_15.20$Shannon ~ Season, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$Shannon ~ Depth, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$Shannon ~ Location, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$Shannon ~ Dataset, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$Shannon ~ Pollution, data = alpha_df_15.20)

kruskal.test(alpha_df_15.20$InvSimpson ~ Season, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$InvSimpson ~ Depth, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$InvSimpson ~ Location, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$InvSimpson ~ Dataset, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$InvSimpson ~ Pollution, data = alpha_df_15.20)

pairwise.wilcox.test(alpha_df_15.20$InvSimpson, alpha_df_15.20$Season,
                     p.adjust.method = "BH")


#####################################
# Explore results on phylum, class and family level
#####################################

#Compare on Phylum level
prevdf_sum_phylum <- plyr::ddply(prevdf3, "Phylum", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))


#Explore results on class level
prevdf_sum_class <- plyr::ddply(prevdf3, "Class", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#Explore results on genus level
prevdf_sum_genus <- plyr::ddply(prevdf3, "Genus", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top class
top_class <- prevdf_sum_class %>% 
  top_n(10, Proportion)
top_class

##############################
#Preparation of data for taxonomic comparison
##############################

#Abundance value transformation
phy_obj3ra = transform_sample_counts(phy_obj3, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra)
write.table(phy_obj3ra.melt, "./tables/phy_obj3ra.melt.15_20.txt")


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


###########################
#Taxonomic compositions on class level
###########################
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
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
phy_obj3_melt.agg.class$Season = factor(phy_obj3_melt.agg.class$Date, levels=c( '15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Class3_Date.pdf", last_plot())


###Subsets by taxonomy - Campylobacteria (position stack/position fill)
Campyl.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Campylobacteria")

Campyl.melt$Date_f = factor(Campyl.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
ggplot(Campyl.melt, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Date_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Campylobacteria.pdf", last_plot())


###Subsets by taxonomy - Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Alphaproteobacteria")
Alpha.melt$Date_f = factor(Alpha.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

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
  facet_grid(Date_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Alphaproteobacteria.pdf", last_plot())

###Subsets by taxonomy - Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Gammaproteobacteria")
Gamma.melt$Date_f = factor(Gamma.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))


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
  facet_grid(Date_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Gammaproteobacteria.pdf", last_plot())


###Subsets by taxonomy - Bacteroidia
Bacter.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacteroidia")
Bacter.melt$Date_f = factor(Bacter.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Bacter.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Bacilli
Bacilli.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Bacilli")
Bacilli.melt$Date_f = factor(Bacilli.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Bacilli.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

###Subsets by taxonomy - Clostridia
Clost.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Clostridia")
Clost.melt$Date_f = factor(Clost.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Clost.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Clostridia.pdf", last_plot())

###Vibrio
Vibrio<-subset(phy_obj3_melt.agg.genus, subset = Genus == "Vibrio")

##################################
#Microbial pollution indicators
##################################

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)
write.table(ps_Mic_Ind_ra.melt, "./tables/ps_Mic_Ind_ra.melt.txt")

#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]

ps_Mic_Ind_ra.melt.agg.genus$Date_f = factor(ps_Mic_Ind_ra.melt.agg.genus$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
colourCount = length(unique(ps_Mic_Ind_ra.melt.agg.genus$Family))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
ggplot(ps_Mic_Ind_ra.melt.agg.genus, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
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
phy_obj3.dds <- phyloseq_to_deseq2(phy_obj3, ~1)
geoMeans = apply(counts(phy_obj3.dds), 1, gm_mean)

# DESeq2 Variance Stabilization
phy_obj3.dds = estimateSizeFactors(phy_obj3.dds, geoMeans = geoMeans)
phy_obj3.dds <- estimateDispersions(phy_obj3.dds)
otu.vst <- getVarianceStabilizedData(phy_obj3.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(phy_obj3))
ps.vst<-phy_obj3
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
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Depth+Dataset, data= df, perm = 999)
adonis_all

#Post-hoc test (??)

df <- as(sample_data(df_ps.ind.vst), "data.frame")
d <- phyloseq::distance(df_ps.ind.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Depth+Dataset, data= df, perm = 999)
adonis_all

#Post-hoc test (??)

#Compare only surface
phy_obj3s <- subset_samples(phy_obj3, Depth == "surface")

# DESeq conversion 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
phy_obj3.dds <- phyloseq_to_deseq2(phy_obj3s, ~1)
geoMeans = apply(counts(phy_obj3.dds), 1, gm_mean)

# DESeq2 Variance Stabilization
phy_obj3.dds = estimateSizeFactors(phy_obj3.dds, geoMeans = geoMeans)
phy_obj3.dds <- estimateDispersions(phy_obj3.dds)
otu.vst <- getVarianceStabilizedData(phy_obj3.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(phy_obj3s))
ps.vst<-phy_obj3s
otu_table(ps.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

# Do the ordination with variance stabilized data
phy_obj3.vst.nmds <- ordinate(ps.vst, "NMDS", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.nmds, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/nMDS_var.pdf")

phy_obj3.vst.rda <- ordinate(ps.vst, "RDA", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.rda, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
ggsave("./output_graphs/RDA_var.pdf")
plot_ordination(ps.vst, phy_obj3.vst.rda, type="taxa", color="Phylum", title="taxa")

df <- as(sample_data(ps.vst), "data.frame")
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Dataset, data= df, perm = 999)
adonis_all


