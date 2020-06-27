##############################
#Phyloseq Microbiome 2015
##############################


#Loading packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")
library(dplyr)
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(vegan); packageVersion("vegan")
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


##################################
#Generate an overview table
##################################

summary_table <- data.frame()

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

ps_chl
ps_mit

pruned_seq_sums <- data.frame(Chloroplast = ps_chl,Mitochondria = ps_mit)
pruned_seq_sums$Row.names <- paste("X",row.names(pruned_seq_sums), sep ="")

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
phy_obj0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") &!Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

#alpha diversity indeces
ps0_alpha <- estimate_richness(phy_obj0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Row.names<-rownames(ps0_alpha)

#merge all together
ps0_summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(nonchim/input,2),
         Chloroplast.Seq.prop = round(Chloroplast/nonchim,2),
         Mitochondria.Seq.prop = round(Mitochondria/nonchim,2)) %>%
  select("Location","Season", #metadata
         "input","filtered", "merged", "nonchim", #dada2 
         "Chloroplast","Mitochondria", #taxa
         "Observed","Chao1","Shannon","InvSimpson") #alpha div

summary_table <- rbind(summary_table,ps0_summary_table)

write.table(summary_table, "./Microbiome_overview_table.txt" , sep = "\t", quote = F)

#Statistical comparison by type of samples

Chao1_Wilcox <- summary_table   %>%
  group_by(Season) %>%
  rstatix::wilcox_test(Chao1 ~ Location, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

Shannon_Wilcox <- summary_table   %>%
  group_by(Season) %>%
  rstatix::wilcox_test(Shannon ~ Location, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

InvSimpson_Wilcox <- summary_table   %>%
  group_by(Season) %>%
  rstatix::wilcox_test(InvSimpson ~ Location, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

#####################################
#Plot rarefaction
####################################



#Alpha diversity
plot_richness(phy_obj0, x="Season", measures=c("Shannon", "Simpson"), color="Location")
ggsave("AlphaDiversity_seasons.pdf", last_plot(), device = "pdf")

plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"), color="Season")
ggsave("AlphaDiversity_location.pdf", last_plot(), device = "pdf")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

p <- plot_ordination(ps.prop, ord.nmds.bray, color="Season", shape="Location", title="Bray NMDS")+ geom_point(size=4)
plot(p)
ggsave("NMDS_Seasons2.pdf", last_plot())



plot_ordination(ps.prop, ord.nmds.bray, color="Season", title="Bray NMDS")
ggsave("NMDS_Location.pdf", last_plot())

#Barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Season", fill="Family") + facet_wrap(~Location, scales="free_x")

top40 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:40]
ps.top40 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top40 <- prune_taxa(top40, ps.top40)
plot_bar(ps.top40, x="Season", fill="Family") + facet_wrap(~Location, scales="free_x")



plot_bar(ps.top30, x="Season", fill="Family") + facet_wrap(~Location, scales="free_x")

#Marge ASV at the class level
#NArm set FALSE forces the function to keep the unclassified OTU at the class level
ps.class = tax_glom(ps, taxrank="Class", NArm=FALSE)
ps.class
plot_bar(ps.class, fill="Class", x="Location")+facet_wrap(~Season, scales="free_x", nrow=1)
ps.order = tax_glom(ps, taxrank="Order", NArm=FALSE)
ps.order
plot_bar(ps.order, fill="Order", x="Location")+facet_wrap(~Season, scales="free_x", nrow=1)
top50.class <- names(sort(taxa_sums(ps.order), decreasing=TRUE))[1:50]
ps.top50.class <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50.class <- prune_taxa(top50.class, ps.top50.class)
plot_bar(ps.top50.class, x="Location", fill="Order") + facet_wrap(~Season, scales="free_x", nrow=1)


#Top 1000
top1000 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:1000]
ps.top1000 <- prune_taxa(top1000, ps)
ps.top1000.class = tax_glom(ps.top1000, taxrank="Class", NArm=FALSE)
ps.top1000.class
plot_heatmap(ps.class, "NMDS", "bray", "Location", "Class", low="#FFFFCC", high="#000033", na.value="white")+facet_wrap(~Season, scales="free_x", nrow=1)


plot_bar(ps.class, fill="Class", x="Location")+facet_wrap(~Season, scales="free_x", nrow=1)
ps.order = tax_glom(ps, taxrank="Order", NArm=FALSE)
ps.order

#HeatMap
gpt <- subset_taxa(ps, Kingdom=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:30]), gpt)
plot_heatmap(gpt, sample.label="Location")

p <- plot_heatmap(gpt, "NMDS", "bray", "Location", "Family", low="#FFFFCC", high="#000033", na.value="white")
p

p.class <- plot_heatmap(ps.top50.class, "NMDS", "bray", "Location", "Class", low="#FFFFCC", high="#000033", na.value="white")
p.class


#OTU that represent at least 20% of reads in at least one sample
ps.top <- filter_taxa(ps, function(x) sum(x > total*0.20) > 0, TRUE)
ps.top



#info
nsamples(ps)
ntaxa(ps)
#total count in each sample
head(sample_sums(ps))
rank_names(ps)

