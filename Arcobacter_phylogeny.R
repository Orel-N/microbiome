##############################################
#Arcobacter phylogenetic tree
##############################################

#Load packages
library(dplyr)
library(ggplot2)
library(phyloseq)
library(Biostrings)
library("ape")
library("ggtree") 
library(tidyr)


#Load data
phy_obj3 <- readRDS("./phyloseqFiltered.RDS")

##Both datasets (2015, 2020)
phy_obj3_2 <- phy_obj3
phy_obj3_2.ra <- transform_sample_counts(phy_obj3_2, function(x) x / sum(x))

##Only 2020
phy_obj3 <- subset_samples(phy_obj3, Dataset == "2020")
phy_obj3 <- prune_taxa(taxa_sums(phy_obj3)>0, phy_obj3)
phy_obj3

phy_obj3.ra <- transform_sample_counts(phy_obj3, function(x) x / sum(x))


##################################
#Subset taxa and extract seqences
##################################

#Arcobacter family

##Only 2020
Arcobacter_family <- subset_taxa(phy_obj3.ra, Family == "Arcobacteraceae")
saveRDS(Arcobacter_family, "./data/Arcobacter_family_ps.RDS")

Arcobacter_family.rs <- refseq(Arcobacter_family)
writeXStringSet(Arcobacter_family.rs, "./Arcobacter_family.fasta", format = "fasta")

##Both datsets
Arcobacter_family_2 <- subset_taxa(phy_obj3_2.ra, Family == "Arcobacteraceae")
saveRDS(Arcobacter_family_2, "./data/Arcobacter_family_ps_15_20.RDS")

Arcobacter_family_2.rs <- refseq(Arcobacter_family_2)
writeXStringSet(Arcobacter_family_2.rs, "./Arcobacter_family_15_20.fasta", format = "fasta")



#Similar for other groups

Campylobacteria_class <- subset_taxa(phy_obj3.ra, Class == "Campylobacteria")
Vibrio_family <- subset_taxa(phy_obj3.ra, Family == "Vibrionaceae")
Mycobacteriaceae_family <- subset_taxa(phy_obj3.ra, Family == "Mycobacteriaceae")

Campylobacteria_class.rs <- refseq(Campylobacteria_class)
writeXStringSet(Campylobacteria_class.rs, "./Campylobacteria_class.fasta", format = "fasta")

Vibrio_family.rs <- refseq(Vibrio_family)
writeXStringSet(Vibrio_family.rs, "./Vibrio_family.fasta", format = "fasta")

Mycobacteriaceae_family.rs <- refseq(Mycobacteriaceae_family)
writeXStringSet(Mycobacteriaceae_family.rs, "./Mycobacteriaceae_family.fasta", format = "fasta")

Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
ps_Mic <- subset_taxa(phy_obj3.ra, Family %in% c(Mic_Ind$Family))


###################################
#Phylogenetic tree with closest relatives (SILVA platform)
###################################

#Upload Arcobacter_family.fasta on Silva webside and create tree 
#(Classify 0.97 similarity, compute tree: Denovo including neighbours, program:RAxML, model: GTR)
#View in wasabi and download results

#Import results in R
filename <- "./data/Arcobacter.nex"
t1 <- ape::read.nexus(filename)
t1
plot_tree(t1, label.tips="taxa_names")

ggtree(t1, layout = "circular") + geom_tiplab()
ggtree(t1) + geom_tiplab(as_ylab=TRUE, hjust = -.2, color='firebrick')

#Problem: 
#Show names of closest neighbours
#Some closest neighbours are "uncultured bacterium", "uncultured Epsilonproteobacteria"..


############################################
#Create phylogenetic tree of phyloseq object
###########################################

###Create tree for Arcobacter phyloseq

random_tree = rtree(ntaxa(phy_obj3.ra), rooted=TRUE, tip.label=taxa_names(phy_obj3.ra))
phy_obj3ra.2 = merge_phyloseq(phy_obj3.ra, random_tree)

Arcobacter_family_2 <- subset_taxa(phy_obj3ra.2, Family == "Arcobacteraceae")
Arcobacter_family_2

plot_tree(Arcobacter_family_2, label.tips="taxa_names", color = "Location",
          size="abundance", ladderize = "left")

ggtree(Arcobacter_family_2) + geom_tiplab() 



##Both datasets
random_tree = rtree(ntaxa(phy_obj3_2.ra), rooted=TRUE, tip.label=taxa_names(phy_obj3_2.ra))
phy_obj3ra_2.2 = merge_phyloseq(phy_obj3_2.ra, random_tree)

class(sample_data(phy_obj3ra_2.2)$Dataset)
sample_data(phy_obj3ra_2.2)$Dataset <- as(sample_data(phy_obj3ra_2.2)$Dataset, "character")

Arcobacter_family_2_15_20 <- subset_taxa(phy_obj3ra_2.2, Family == "Arcobacteraceae")
Arcobacter_family_2_15_20


plot_tree(Arcobacter_family_2_15_20, label.tips="taxa_names", color = "Location", shape = "Dataset",
          size="abundance", ladderize = "left")


sample_data(Dataset)
ggtree(Arcobacter_family_2_15_20) + geom_tiplab() 



###Create phylogenetic tree for Arcobacter phyloseq with reference seq

##Create phyloseq with reference seq

#Find reference seq for Arcobacter in studies: (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4224126/ 
#and https://www.sciencedirect.com/science/article/pii/S0043135419311790 )
#Download FASTA from NCBI
#Help for the code: https://github.com/joey711/phyloseq/issues/1150 


#Load FASTA of ref seq in R
rs <- readDNAStringSet(file="./data/Arcobacter_ref_seq_NCBI.fasta", format = "fasta", 
                       nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
seq_name = names(rs)
sequence = paste(rs)
df <- data.frame(seq_name, sequence)
df

#Create otu_table
otumat <- matrix(0.0001, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

#Create metadata table
SAM <- sample_data(Arcobacter_family)[1,]
SAM[,] <- NA
SAM$Dataset <- "Reference"
sample_names(SAM) <- "Reference"
# or : rownames(SAM) <- "Reference"

#Create taxonomy table

tax_df <- data.frame(wholetax = rep(c("level1.level2.level3.level4.level5.level6"), 36))
rownames(tax_df) <- names(rs)  
tax_df <- tax_df %>% separate(wholetax, c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6"))
ref_taxa <- tax_table(as.matrix(tax_df))
ref_taxa

#Create phyloseq
ps.ref <- phyloseq(OTU, SAM, rs, ref_taxa)
ps.ref

##Merge Arcobacter phyloseq and reference phyloseq: Only 2020

ps.merged <- merge_phyloseq(Arcobacter_family, ps.ref)

ps.merged
Arcobacter_family
ps.ref

##Calculate phylogenetic tree

#load libraries
library(DECIPHER)
library(phangorn)

# calculate tree
seqs <- refseq(ps.merged)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
phang.align <- as.phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload = TRUE)

# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(fitGTR$tree))
ps.merged

##Plot tree
plot_tree(ps.merged, label.tips="taxa_names", ladderize = "left")

plot_tree(ps.merged, label.tips="taxa_names", color = "Location", size = "abundance", ladderize = "left")


##Merge Arcobacter phyloseq and reference phyloseq: Both datasets

ps.merged <- merge_phyloseq(Arcobacter_family_2, ps.ref)

ps.merged
Arcobacter_family_2
ps.ref

##Calculate phylogenetic tree
#https://f1000research.com/articles/5-1492/v2 

#load libraries
library(DECIPHER)
library(phangorn)

# calculate tree
seqs <- refseq(ps.merged)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
phang.align <- as.phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload = TRUE)

# merge ps object with tree
ps.merged <- merge_phyloseq(ps.merged, phy_tree(fitGTR$tree))
ps.merged

##Plot tree
plot_tree(ps.merged, label.tips="taxa_names", ladderize = "left")

plot_tree(ps.merged, label.tips="taxa_names", text.size = 2.5, color = "Location", size = "abundance", shape = "Dataset", base.spacing=0.03, ladderize = "left")
