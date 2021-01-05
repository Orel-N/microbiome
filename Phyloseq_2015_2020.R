##############################
#Phyloseq Microbiome 2015_2020
##############################

#Load packages
library(dplyr)
library(ggplot2)
library(phyloseq)

#Set theme
theme_set(theme_bw())


##################################
#Import dada2 output into phyloseq
##################################

#Import "ASV table" (matrix)
seqtab.nochim <- read.csv("./dada2/dada2_seqtab_nochim2_2.txt", h=T, sep="\t")

#Import taxonomy table (matrix)
taxa <- as.matrix(read.csv("./dada2/dada2_taxonomy_table_2.txt", h=T,sep = "\t"))

#Import Sample Data - "metadata"
meta <- read.csv("./data/Metadata_2015_2020.csv", h=T, row.names=1, sep = ",")

#Edit metadata
rownames(meta)<- paste("X",meta$SampleID, sep ="")

#Check order of seq in taxonomy table and in seq table
all.equal(rownames(seqtab.nochim), rownames(taxa))

#Create a phyloseq object from the ASV table and taxonomy assigned by DADA2
ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), 
              sample_data(meta), tax_table(taxa))
ps

#add reference sequence and replace variants with ASVs
ps2 <- Biostrings::DNAStringSet(taxa_names(ps))
names(ps2) <- taxa_names(ps)

ps <- merge_phyloseq(ps, ps2)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


####################
#Taxonomic filtering
####################

#Show available ranks in the dataset
rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

#check how many ASVs were unclassified on phylum level, or assigned to Eukaryota, Chloroplast, Mitochondria and Archaea
ps_euk <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Eukaryota"))),"vector")
ps_chl <- as(sample_sums(subset_taxa(ps, Order %in% c("Chloroplast"))),"vector")
ps_mit <- as(sample_sums(subset_taxa(ps, Family %in% c("Mitochondria"))),"vector")
ps_arch <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Archaea"))),"vector")
ps_all <- as(sample_sums(subset_taxa(ps)),"vector")

#remove unclassified on phylum level, chloroplast, Mitochondrial and Archaeal sequence variants
ps0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") & !is.na(Phylum) & !Phylum %in% c("", "uncharacterized") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") & !Kingdom %in% c("Archaea"))
ps0
#Create separate phyloseq object with seq assigned to Chloroplast
ps0_chl <- subset_taxa(ps, Order == "Chloroplast")
ps0_chl

#check how many ASVs remains after filtering
ps_all_after_removal <- as(sample_sums(subset_taxa(ps0)),"vector")


######################
#Create overview table
######################

#Add Row.names column in metadata
meta$Row.names<-paste(row.names(meta))

#Summary of dada2 workflow
reads.tab <- read.csv2(file.path("./dada2/Overview_dada2_15_20_2.txt"),header = TRUE, sep = "\t")
reads.tab$Row.names <-paste("X",row.names(reads.tab), sep ="")

#Create overview table of taxonomic filtration
pruned_seq_sums <- data.frame(Eukaryota = ps_euk, Chloroplast = ps_chl, Mitochondria = ps_mit, Archaea = ps_arch, All_phyloseq = ps_all, After_tax_filt = ps_all_after_removal)
pruned_seq_sums$Row.names <- paste(row.names(pruned_seq_sums))
pruned_seq_sums

#Alpha diversity calculation
ps0_alpha <- estimate_richness(ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Evenness <- ps0_alpha$Shannon/log(ps0_alpha$Observed)
ps0_alpha$Row.names<-rownames(ps0_alpha)

#Merge together overview data: metadata, dada2 overview, taxonomy filtration overview and alpha diversity 
summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(nochim/input,2),
         Euk = round(Eukaryota,2),
         Chloroplast = round(Chloroplast,2),
         Mitochondria = round(Mitochondria,2),
         Archaea = round(Archaea,2)) %>%
  select("SampleID","Location","Season", "Depth", "Dataset", #metadata
         "input","filtered", "merged", "nochim", "Seq.prop", #dada2 overview
         "Chloroplast","Mitochondria", "Archaea", "After_tax_filt",#taxonomy filtration overview
         "Observed","Chao1","Shannon","InvSimpson", "Evenness") #alpha diversity

write.table(summary_table, "./output_tables/16S_15_20_overview_table.txt" , sep = "\t", quote = F)


############################################
#ASVs distribution and prevalence filtering
############################################

#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
prev = apply(X = otu_table(ps0),
             MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

#Plot prevalence
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts_BeforeF.pdf", last_plot())

#Create a table to compare on phylum level
prev_sum_phylum <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion=100*Abundance/sum(Abundance))
prev_sum_phylum

write.table(prev_sum_phylum, "./output_tables/prev_sum_phylum.txt", sep = "\t", quote = F)

#Define phyla to filter (abundance on phylum level less than 100)
fliterPhyla = subset(prev_sum_phylum, Abundance < 100)

#filter phyla and save phyloseq
ps3 = subset_taxa (ps0, !Phylum %in% fliterPhyla$Phylum)
ps3

#Plot prevalence
prevdf3 = subset(prevdf, Phylum %in% get_taxa_unique(ps3, "Phylum"))
ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(ps3),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") 

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts_AfterF.pdf", last_plot())

#Save phyloseq objects
saveRDS(ps, "./data/phyloseqRaw.RDS")
saveRDS(ps0, "./data/phyloseqTaxFiltered.RDS")
saveRDS(ps3, "./data/phyloseqPrevFiltered.RDS")


# CLEAN UP #################################################

# Clear plots
dev.off()

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
