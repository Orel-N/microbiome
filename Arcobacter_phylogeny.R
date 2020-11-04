
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(RColorBrewer)
library(iNEXT)
library(DESeq2)
library(rstatix)
library(vegan)
library(Biostrings)


phy_obj3 <- readRDS("./phyloseqFiltered.RDS")
prevdf3 <- readRDS("./prevdfFiltered.RDS")

##############################
#Preparation of data for taxonomic comparison
##############################

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


#Calculate abundance for each Class
phy_obj3_melt.agg.class <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.class$Abundance <- phy_obj3_melt.agg.class$Abundance*100
phy_obj3_melt.agg.class<- phy_obj3_melt.agg.class[phy_obj3_melt.agg.class$Abundance>0,]



Campylobacteria_class <- subset_taxa(phy_obj3, Class == "Campylobacteria")
Arcobacter_family <- subset_taxa(phy_obj3, Family == "Arcobacteraceae")
Vibrio_family <- subset_taxa(phy_obj3, Family == "Vibrionaceae")
Mycobacteriaceae_family <- subset_taxa(phy_obj3, Family == "Mycobacteriaceae")

Campylobacteria_class.rs <- refseq(Campylobacteria_class)
writeXStringSet(Campylobacteria_class.rs, "./Campylobacteria_class.fasta", format = "fasta")

Arcobacter_family.rs <- refseq(Arcobacter_family)
writeXStringSet(Arcobacter_family.rs, "./Arcobacter_family.fasta", format = "fasta")

Vibrio_family.rs <- refseq(Vibrio_family)
writeXStringSet(Vibrio_family.rs, "./Vibrio_family.fasta", format = "fasta")

Mycobacteriaceae_family.rs <- refseq(Mycobacteriaceae_family)
writeXStringSet(Mycobacteriaceae_family.rs, "./Mycobacteriaceae_family.fasta", format = "fasta")

##################################
#Microbial pollution indicators
##################################

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic_Ind_ra <- subset_taxa(phy_obj3ra, Family %in% c(Mic_Ind$Family))

saveRDS(ps_Mic_Ind_ra, "./phyloseq_Mic_Ind_ra")

library(ape)
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree") 
library(phyloseq)
BiocManager::install("ggtree")
t <- read.tree("./data/Arcobacter.tree")
t

ggtree(t) + geom_text2(t, )
plot.phylo(t)
plot(t1, type = "phylogram") 
p

filename <- "./data/Arcobacter.nex"
t1 <- ape::read.nexus(filename)
t1
plot_tree(t1, label.tips="taxa_names")
ggtree(t1, layout = "circular") + geom_tiplab()
ggtree(t1) + geom_tiplab(as_ylab=TRUE, hjust = -.2, color='firebrick')

#Tree
library("ape")
library(ggjoy)
install.packages('ggjoy')

random_tree = rtree(ntaxa(phy_obj3ra), rooted=TRUE, tip.label=taxa_names(phy_obj3ra))
plot(random_tree)
phy_obj3ra.2 = merge_phyloseq(phy_obj3ra, random_tree)

Arcobacter_family_2 <- subset_taxa(phy_obj3ra.2, Family == "Arcobacteraceae")
plot_tree(Arcobacter_family_2, label.tips="taxa_names", color = "Location",
          size="abundance", ladderize = "left")


Arcobacter_family_2

ggtree(Arcobacter_family_2) + geom_tiplab() 

