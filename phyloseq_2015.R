##############################
#Phyloseq Microbiome 2015
##############################


#Loading packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(vegan); packageVersion("vegan")
library(iNEXT); packageVersion("iNEXT")

theme_set(theme_bw())

##################################
#Import dada2 output into phyloseq
##################################

#Import Sample Data - "metadata"(data.frame)
meta <- read.table("./data/metadata.txt", row.names = 1)
meta

#Import "ASV table" (matrix)
seqtab.nochim <- readRDS ("./data/seqtab.nochim.rds")
#correct names
#colnames(seqtab.nochim) <- paste0("ASV", 1:ncol(seqtab.nochim))

#Import taxonomy table (matrix)
taxa <- readRDS ("./data/taxa.rds")

#Check order
all.equal(colnames(seqtab.nochim), rownames(taxa))

#Create a phyloseq object from the OTU table/ASV table and taxonomy assigned by DADA2
ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa))
ps

#add reference sequence and replace variants with ASVs
ps2 <- Biostrings::DNAStringSet(taxa_names(ps))
names(ps2) <- taxa_names(ps)
ps <- merge_phyloseq(ps, ps2)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#remove unobserved ASVs
ps<- prune_taxa(taxa_sums(ps)>0, ps)
ps

##################################
#Generate an overview table
##################################

#metadata
meta <- as(sample_data(ps), "data.frame")
meta$Row.names<-paste("X",row.names(meta), sep ="")

#summary of dada2 workflow
track <- readRDS ("./data/track.rds")
reads.tab <- as.data.frame(track)
reads.tab$Row.names <-paste("X",row.names(reads.tab), sep ="")

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
ps_euk <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Eukaryota"))),"vector")
ps_uncl <- as(sample_sums(subset_taxa(ps, Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl"))),"vector")
ps_chl <- as(sample_sums(subset_taxa(ps, Order %in% c("Chloroplast"))),"vector")
ps_mit <- as(sample_sums(subset_taxa(ps, Family %in% c("Mitochondria"))),"vector")
ps_euk <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Eukaryota"))),"vector")
ps_arch <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Archaea"))),"vector")

ps_chl
ps_mit
ps_arch

pruned_seq_sums <- data.frame(Chloroplast = ps_chl,Mitochondria = ps_mit, Archaea = ps_arch)
pruned_seq_sums$Row.names <- paste("X",row.names(pruned_seq_sums), sep ="")

#remove unclassified on phylum level, chloroplast, Mitochondrial and Archaeal sequence variants
phy_obj0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") &!Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") & !Kingdom %in% c("Archaea"))

#alpha diversity indeces
ps0_alpha <- estimate_richness(phy_obj0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Row.names<-rownames(ps0_alpha)

#merge all together
summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  select("Location","Season", #metadata
         "input","filtered", "merged", "nonchim", #dada2 
         "Chloroplast","Mitochondria", "Archaea", #taxa
         "Observed","Chao1","Shannon","InvSimpson") #alpha div


#save the summary table

write.table(summary_table, "./Micro_overview_table.txt" , sep = "\t", quote = F)



#####################################
#Plot rarefaction
####################################

phy_obj <- phy_obj0

#https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

ps0.iNEXT <- iNEXT(as.data.frame(otu_table(phy_obj)), q=0, datatype="abundance", knots = 40)

ggiNEXT(ps0.iNEXT, type = 1)

