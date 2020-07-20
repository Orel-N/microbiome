##############################
#Phyloseq Microbiome 2015
##############################
#instal phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('phyloseq')

install.packages("dplyr")
install.packages("rstatix")
install.packages("iNEXT")

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
reads.tab

#Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
ps_euk <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Eukaryota"))),"vector")
ps_uncl <- as(sample_sums(subset_taxa(ps, Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl"))),"vector")
ps_chl <- as(sample_sums(subset_taxa(ps, Order %in% c("Chloroplast"))),"vector")
ps_mit <- as(sample_sums(subset_taxa(ps, Family %in% c("Mitochondria"))),"vector")
ps_arch <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Archaea"))),"vector")


ps_chl
ps_mit
ps_arch

pruned_seq_sums <- data.frame(Chloroplast = ps_chl,Mitochondria = ps_mit, Archaea = ps_arch)
pruned_seq_sums$Row.names <- paste("X",row.names(pruned_seq_sums), sep ="")
pruned_seq_sums

#remove unclassified on phylum level, chloroplast, Mitochondrial and Archaeal sequence variants
ps0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") &!Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") & !Kingdom %in% c("Archaea"))


#alpha diversity indeces
ps0_alpha <- estimate_richness(ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Row.names<-rownames(ps0_alpha)

#Edit number of decimals
ps0_alpha$Chao1 <- round(ps0_alpha$Chao1, 2)
ps0_alpha$se.chao1 <- round(ps0_alpha$se.chao1, 2)
ps0_alpha$Shannon <- round(ps0_alpha$Shannon, 2)
ps0_alpha$InvSimpson <- round(ps0_alpha$InvSimpson, 2)
ps0_alpha



#merge all together
summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  select("Location","Season", #metadata
         "input","filtered", "merged", "nonchim", #dada2 
         "Chloroplast","Mitochondria", "Archaea", #taxa
         "Observed","Chao1","Shannon","InvSimpson") #alpha div


#save the summary table

write.table(summary_table, "./output_tables/MT_overview_table.txt" , sep = "\t", quote = F)



####################################
#Plot rarefaction
####################################

#https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

#transpose phyloseq object (speciesXsample -> sampleXspecies)
phy_obj <- t(ps0)
phy_obj

ps0.iNEXT <- iNEXT(as.data.frame(otu_table(phy_obj)), q=0, datatype="abundance", knots = 40)

ggiNEXT(ps0.iNEXT, type = 1)
ggsave("./output_graphs/MT_rarefactions.pdf", plot=last_plot())

ggiNEXT(ps0.iNEXT, type = 2)
ggsave("./output_graphs/MT_rarefactions_coverage.pdf", plot=last_plot())


##################################
#Filtering
###################################

#Taxonomic filtering
####################

# Create table, number of features for each phyla
t2 <- table(tax_table(phy_obj)[, "Phylum"], exclude = NULL)
write.table(t2, "NoFeaturesPhyla.csv")

#Remove features with a Phylum NA
phy_obj
phy_obj2 <- subset_taxa(phy_obj, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
phy_obj2

plot_bar(phy_obj2, fill ="Phylum")
plot_bar(phy_obj, fill = "Phylum ")

#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
# Compute prevalence of each feature, store as data.frame
prev = apply(X = otu_table(phy_obj2),
             MARGIN = ifelse(taxa_are_rows(phy_obj2), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(phy_obj2),
                    tax_table(phy_obj2))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


#Define phyla to filter (observed only one feature)

fliterPhyla = c("Abditibacteriota", "Latescibacterota", "PAUC34f", "Schekmanbacteria", "NB1-j", "Nitrospirota", "Hydrogenedentes")
#filter phyla
phy_obj3 = subset_taxa (phy_obj2, !Phylum %in% fliterPhyla)
plot_bar(phy_obj3, fill ="Phylum")

#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
prevdf3 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj3, "Phylum"))



#Prevalence filtering
########################
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(phy_obj2),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") 

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts.pdf", last_plot())

# Subset to the remaining phyla
prevdf3 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj3, "Phylum"))
ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(phy_obj3),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts2.pdf", last_plot())


#remove singletons
prevdf1 <- prevdf %>% filter (TotalAbundance>1)
prevdf1
prevdf2 <- prevdf %>% filter (Prevalence>1)
prevdf2

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phy_obj2),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ggplot(prevdf2, aes(TotalAbundance, Prevalence / nsamples(phy_obj2),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


#####################################
# Explore results on phylum, class and family level
#####################################

#Compare on Phylum level
prevdf_sum_phylum <- plyr::ddply(prevdf, "Phylum", function(prevdf){cbind(Samples=mean(prevdf$Prevalence),Abundance=sum(prevdf$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
prevdf_sum_phylum3 <- plyr::ddply(prevdf3, "Phylum", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top phyla
top_Phylum <- prevdf_sum_phylum %>% 
  top_n(10, Proportion)
top_Phylum

#Explore results on class level
prevdf_sum_class <- plyr::ddply(prevdf3, "Class", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
#define top class
top_class <- prevdf_sum_class %>% 
  top_n(10, Proportion)
top_class


prevdf_sum_family <- plyr::ddply(prevdf3, "Family", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
#define top class
top_family <- prevdf_sum_family %>% 
  top_n(50, Proportion)
top_family


##############################
#Preparation of data for taxonomic comparison
##############################
#Abundance value transformation
phy_obj3ra = transform_sample_counts(phy_obj3, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra)

#Calculate abundance for each taxa
phy_obj3_Genus <- as.data.frame(as.list(aggregate(Abundance~Location+Season+Class+Order+Family+Genus, phy_obj3ra.melt,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_Genus$Abundance.sum <- phy_obj3_Genus$Abundance.sum*100
phy_obj3_Genus<- phy_obj3_Genus[phy_obj3_Genus$Abundance.sum>0,]


###########################
#Taxonomic compositions on class level
###########################
#Calculate abundance for each Class
phy_obj3_melt.agg.class <- aggregate (Abundance~Location+Season+Class, phy_obj3ra.melt, FUN="sum")
phy_obj3_melt.agg.class$Abundance <- phy_obj3_melt.agg.class$Abundance*100
phy_obj3_melt.agg.class<- phy_obj3_melt.agg.class[phy_obj3_melt.agg.class$Abundance>0,]

#remove below 1% ra
threshold<- 1
phy_obj3_melt.agg.class$Class <- as.character(phy_obj3_melt.agg.class$Class)

taxa_classes <- unique(phy_obj3_melt.agg.class$Class[!phy_obj3_melt.agg.class$Abundance<threshold])

phy_obj3_melt.agg.class$Class[phy_obj3_melt.agg.class$Abundance<threshold] <- "Other taxa"

phy_obj3_melt.agg.class$Class <- factor(phy_obj3_melt.agg.class$Class,
                                         levels=c(taxa_classes,"Other taxa"))

#Plot Class(Location)
phy_obj3_melt.agg.class$Location_f = factor(phy_obj3_melt.agg.class$Location, levels=c('00RI','ERI2','000K'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Season, fill = Class))+
  facet_grid(Location_f~., space= "fixed")+
  geom_bar(stat = "identity", position="fill") +
  scale_y_discrete(limits = c('autumn', 'summer','spring','winter'))
ggsave("./output_graphs/Class2_location.pdf", last_plot())


#Plot Class(Season)
phy_obj3_melt.agg.class$Season_f = factor(phy_obj3_melt.agg.class$Season, levels=c('autumn', 'summer','spring','winter'))
ggplot(phy_obj3_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season_f~., space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_y_discrete(limits=c('00RI','ERI2','000K'))
ggsave("./output_graphs/Class2_season.pdf", last_plot())


#calculate abundance for each Order
phy_obj3_order <- as.data.frame(as.list(aggregate(Abundance~Location+Season+Class+Order, phy_obj3ra.melt,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_order$Abundance.sum <- phy_obj3_order$Abundance.sum*100
phy_obj3_order<- phy_obj3_order[phy_obj3_order$Abundance.sum>0,]



#calculate abundance for each Family
phy_obj3_family <- as.data.frame(as.list(aggregate(Abundance~Location+Season+Class+Order+Family, phy_obj3ra.melt,
                                                   FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj3_family$Abundance.sum <- phy_obj3_family$Abundance.sum*100
phy_obj3_family<- phy_obj3_family[phy_obj3_family$Abundance.sum>0,]

#Plot nMDS
phy_obj3.nmds <- ordinate(phy_obj3ra, method = "NMDS", 
                             distance = "jsd")
stressplot(phy_obj3.nmds)

plot_ordination(phy_obj3ra, phy_obj3.nmds, color = "Season")+
  stat_ellipse(type="norm")+
  geom_point(aes(shape = Location),color ="black", size = 5)+
  geom_point(aes(shape = Location), size = 4)+
  coord_fixed()

ggsave("./output_graphs/nMDS.pdf", last_plot())

#Subsets by taxonomy

plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Campilobacterota"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Location",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
psOrd = subset_taxa(phy_obj3ra, Order == "Campylobacterales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL) 

###############
#SESeq2
###############

BiocManager::install("DESeq2")
#Load DESeq2
library("DESeq2")
packageVersion("DESeq2")

#DESeq2 conversion and call
#phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, 
#using the experimental design formula, also shown (the ~Location term)
#The DESeq function does the rest of the testing, in this case with default testing framework
#Comparison on Location 00RI and 000K (ERI2 was excluded)

head(sample_data(phy_obj3)$Location)

phy_obj3a = subset_samples(phy_obj3, Location != "ERI2")

diagdds = phyloseq_to_deseq2(phy_obj3a, ~ Location)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#Investigate test results table
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy_obj3a)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#Plot summary of the results=ASVs that were significantly different between the two Locations
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

##Count how many NAs increase and decrease in abundance
sigtab.pos = sigtab[(sigtab$log2FoldChange >=0),]
pos.row.has.na <- apply(sigtab.pos, 1, function(x){any(is.na(x))})
sum(pos.row.has.na)

sigtab.neg = sigtab[(sigtab$log2FoldChange <=0),]
neg.row.has.na <- apply(sigtab.neg, 1, function(x){any(is.na(x))})
sum(neg.row.has.na)

ggsave("./output_graphs/DESeq2.pdf", last_plot())
