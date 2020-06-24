##############################
#Phyloseq Microbiome 2015
##############################


#Loading packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

###Import dada2 output into phyloseq
#metadata
meta <- read.table("C:/Users/nezao/Documents/1-MIKROBIOM/1-Eksperimentalno delo/2. Eksperiment - TEREN/Rezultati/BacterialCommunity/2015/metadata.txt", row.names = 1)
meta


#Create a phyloseq object from the ASV table and taxonomyassigned by DADA2
ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa))
ps

row.names(seqtab.nochim)
sample_names(seqtab.nochim)


#Alpha diversity
plot_richness(ps, x="Season", measures=c("Shannon", "Simpson"), color="Location")
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

