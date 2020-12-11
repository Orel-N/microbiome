#MICROBIOME - Diversity analysis 2015 and 2020 dataset
#Neža Orel, neza.orel@nib.si
#Script for diversity analyses: graphs and statistic

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
prevdf3 <- readRDS("./prevdfFiltered.RDS")

#Import collor pallets
phyla.col <- readRDS("./data/phyla_col.RDS")
indicators.col <- readRDS("./data/indicators_col.RDS")
season.col <- readRDS("./data/season_col.RDS")
location.col <- readRDS("./data/location_col.RDS")


###################################################
# Explore results on phylum, class and family level
###################################################

#Compare on Phylum level
prevdf_sum_phylum <- plyr::ddply(prevdf3, "Phylum", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#Explore results on class level
prevdf_sum_class <- plyr::ddply(prevdf3, "Class", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top class
top_class <- prevdf_sum_class %>% 
  top_n(10, Proportion)
top_class

#Explore results on genus level
prevdf_sum_genus <- plyr::ddply(prevdf3, "Genus", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))


#############################################
#Preparation of data for taxonomic comparison
#############################################

#Abundance value transformation
phy_obj3ra = transform_sample_counts(phy_obj3, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra)

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


######################################
#Taxonomic compositions on class level
######################################

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
write.csv(phy_obj3_melt.agg.class.table, "./output_tables/Class.csv")

#Plot Class(Date)
phy_obj3_melt.agg.class$Date = factor(phy_obj3_melt.agg.class$Date, levels=c( '4.12.2018','25.04.2019','18.06.2019','15.10.2019','4.02.2016','24.03.2015','15.07.2015','25.11.2015'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Date~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=phyla.col)+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Class3_Date.pdf", last_plot())


###Subsets by taxonomy - Campylobacteria (position stack/position fill)
Campyl.melt <- subset(phy_obj3_melt.agg.family , subset = Class == "Campylobacteria")
ggplot(Campyl.melt, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Date~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))

###Subsets by taxonomy - Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Alphaproteobacteria")

#remove below 1% ra
threshold<- 1
Alpha.melt$Order <- as.character(Alpha.melt$Order)
taxa_classes <- unique(Alpha.melt$Order[!Alpha.melt$Abundance<threshold])
Alpha.melt$Order[Alpha.melt$Abundance<threshold] <- "Other taxa"
Alpha.melt$Order <- factor(Alpha.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

ggplot(Alpha.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))

###Subsets by taxonomy - Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj3_melt.agg.family , subset = Class == "Gammaproteobacteria")

#remove below 1% ra
threshold<- 1
Gamma.melt$Order <- as.character(Gamma.melt$Order)
taxa_classes <- unique(Gamma.melt$Order[!Gamma.melt$Abundance<threshold])
Gamma.melt$Order[Gamma.melt$Abundance<threshold] <- "Other taxa"
Gamma.melt$Order <- factor(Gamma.melt$Order,
                           levels=c(taxa_classes,"Other taxa"))

ggplot(Gamma.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Gammaproteobacteria.pdf", last_plot())


##################################
#Microbial pollution indicators
##################################

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra, Family %in% c(Mic_Ind$Family))

#Melt data
ps_Mic_Ind_ra.melt <- psmelt(ps_Mic_Ind_ra)

#Calculate abundance for each taxa
ps_Mic_Ind_ra.melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, ps_Mic_Ind_ra.melt,
                                                                FUN = function(x) c(sum = sum(x), count=length(x)))))
ps_Mic_Ind_ra.melt.agg.genus$Abundance <- ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum*100
ps_Mic_Ind_ra.melt.agg.genus<- ps_Mic_Ind_ra.melt.agg.genus[ps_Mic_Ind_ra.melt.agg.genus$Abundance.sum>0,]

ps_Mic_Ind_ra.melt.agg.genus$Date_f = factor(ps_Mic_Ind_ra.melt.agg.genus$Date, levels=c( '4.12.2018','25.04.2019','18.06.2019','15.10.2019','4.02.2016','24.03.2015','15.07.2015','25.11.2015'))

ggplot(ps_Mic_Ind_ra.melt.agg.genus, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Date_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_fill_manual(values=indicators.col)+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/MicrobialInidicators.pdf", last_plot())

# Heatmap
ps_Mic_Ind_ra_glom = tax_glom(ps_Mic_Ind_ra, taxrank="Family")
fig <- plot_heatmap(ps_Mic_Ind_ra_glom, "Family", taxa.label = "Family", sample.label="Location", taxa.order="Class", sample.order="Location")
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
phy_obj3.vst.rda <- ordinate(ps.vst, "RDA", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.rda, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
plot_ordination(ps.vst, phy_obj3.vst.rda, label = "Dataset", shape = "Depth", color = "Season") + 
  theme_bw() + geom_point(size=3) + scale_color_manual(values=season.col)

ggsave("./output_graphs/RDA_diversity.pdf")
saveRDS(ps.vst, "./data/ps.vst.RDS")


###########
#PERMANOVA
###########

#statistical significance of the groups
df <- as(sample_data(ps.vst), "data.frame")
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Depth+Dataset, data= df, perm = 999)
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


#Compare only surface
#####################

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
phy_obj3.vst.rda <- ordinate(ps.vst, "RDA", "euclidean")
plot_ordination(ps.vst, phy_obj3.vst.rda, shape = "Location", color = "Season")+facet_wrap(~Depth)+geom_point(size=4)
plot_ordination(ps.vst, phy_obj3.vst.rda, label = "Dataset", color = "Season") + 
  theme_bw() + geom_point(size=3) + scale_color_manual(values=season.col)

#Permanova
df <- as(sample_data(ps.vst), "data.frame")
d <- phyloseq::distance(ps.vst, "euclidean")
adonis_all <- adonis2(d ~ Location+Season+Dataset, data= df, perm = 999)
adonis_all
