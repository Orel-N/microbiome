
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
library("ape")
library("ggtree") 
library(phyloseq)
library(ggjoy)


phy_obj3 <- readRDS("./phyloseqFiltered.RDS")

phy_obj3 <- subset_samples(phy_obj3, Dataset == "2020")
phy_obj3 <- prune_taxa(taxa_sums(phy_obj3)>0, phy_obj3)
phy_obj3

phy_obj3.ra <- transform_sample_counts(phy_obj3, function(x) x / sum(x))


##############################
#Subset taxa
##############################

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

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic <- subset_taxa(phy_obj3, Family %in% c(Mic_Ind$Family))
ps.ra_Mic <- subset_taxa(phy_obj3ra, Family %in% c(Mic_Ind$Family))


##############################
#Silva tree with closest relatives
##############################

filename <- "./data/Arcobacter.nex"
t1 <- ape::read.nexus(filename)
t1
plot_tree(t1, label.tips="taxa_names")
ggtree(t1, layout = "circular") + geom_tiplab()
ggtree(t1) + geom_tiplab(as_ylab=TRUE, hjust = -.2, color='firebrick')


##############################
#Create phylogenetic tree of phyloseq object
##############################


random_tree = rtree(ntaxa(phy_obj3ra), rooted=TRUE, tip.label=taxa_names(phy_obj3ra))
phy_obj3ra.2 = merge_phyloseq(phy_obj3ra, random_tree)

Arcobacter_family_2 <- subset_taxa(phy_obj3ra.2, Family == "Arcobacteraceae")
plot_tree(Arcobacter_family_2, label.tips="taxa_names", color = "Location",
          size="abundance", ladderize = "left")


Arcobacter_family_2

ggtree(Arcobacter_family_2) + geom_tiplab() 

#Create phylogenetic tree with reference seq
library(ape)
library(Biostrings)

rs <- readDNAStringSet("./data/Arcobacter_ref_seq_NCBI.fasta", format = "fasta")
seq_name = names(rs)
sequence = paste(rs)
df <- data.frame(seq_name, sequence)
df

otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

SAM <- sample_data(Arcobacter_family_2)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"
ps.ref <- phyloseq(OTU, SAM, rs)
ps.ref

ps.merged <- merge_phyloseq(Arcobacter_family_2, ps.ref)

ps.merged
Arcobacter_family_2
ps.ref

