##############################
#Phyloseq Microbiome 2015_2020
##############################
#instal phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('phyloseq')

install.packages("dplyr")
install.packages("rstatix")
install.packages("iNEXT")
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(RColorBrewer)

theme_set(theme_bw())

##################################
#Import dada2 output into phyloseq
##################################

#Import Sample Data - "metadata"(data.frame)
meta <- read.table("./data/Metadata_16S_2015_2020_2.txt", header=TRUE)
rownames(meta)<- paste("X",meta$SampleID, sep ="")
meta

#Import "ASV table" (matrix)
seqtab.nochim <- read.csv("./dada2/dada2_seqtab_nochim2_2.txt", h=T, sep="\t")

#Import taxonomy table (matrix)
taxa <- as.matrix(read.csv("./dada2/dada2_taxonomy_table_2.txt", h=T,sep = "\t"))

#Check order
all.equal(rownames(seqtab.nochim), rownames(taxa))
all.equal(colnames(seqtab.nochim), rownames(meta))

#Create a phyloseq object from the OTU table/ASV table and taxonomy assigned by DADA2
ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), sample_data(meta), tax_table(taxa))
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
meta2 <- meta
meta2$Row.names<-paste(row.names(meta))

#summary of dada2 workflow
reads.tab <- read.csv2(file.path("./dada2/Overview_dada2_15_20_2.txt"),header = TRUE, sep = "\t")
reads.tab$Row.names <-paste("X",row.names(reads.tab), sep ="")

#Show available ranks in the dataset
rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
ps_euk <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Eukaryota"))),"vector")
ps_chl <- as(sample_sums(subset_taxa(ps, Order %in% c("Chloroplast"))),"vector")
ps_mit <- as(sample_sums(subset_taxa(ps, Family %in% c("Mitochondria"))),"vector")
ps_arch <- as(sample_sums(subset_taxa(ps, Kingdom %in% c("Archaea"))),"vector")

ps_euk
ps_chl
ps_mit
ps_arch

pruned_seq_sums <- data.frame(Eukaryota = ps_euk,Chloroplast = ps_chl,Mitochondria = ps_mit, Archaea = ps_arch)
pruned_seq_sums$Row.names <- paste("X",row.names(pruned_seq_sums), sep ="")
pruned_seq_sums

#remove unclassified on phylum level, chloroplast, Mitochondrial and Archaeal sequence variants
ps0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") &!Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") & !Kingdom %in% c("Archaea"))
ps0


#alpha diversity indeces
ps0_alpha <- estimate_richness(ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Row.names<-rownames(ps0_alpha)


#merge all together
summary_table <- merge(meta2,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(nochim/input,2),
         Euk.Seq.prop = round(Eukaryota/nochim,2),
         Chloroplast.Seq.prop = round(Chloroplast/nochim,2),
         Mitochondria.Seq.prop = round(Mitochondria/nochim,2),
         Archaea.Seq.prop = round(Archaea/nochim,2)) %>%
  select("Location","Season", "Depth", #metadata
         "input","filtered", "merged", "nochim", "Seq.prop", #dada2 
         "Chloroplast.Seq.prop","Mitochondria.Seq.prop", "Archaea.Seq.prop", #taxa
         "Observed","Chao1","Shannon","InvSimpson") #alpha div

summary_table2 <- merge(meta2,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(nochim/input,2),
         Euk = round(Eukaryota,2),
         Chloroplast = round(Chloroplast,2),
         Mitochondria = round(Mitochondria,2),
         Archaea = round(Archaea,2)) %>%
  select("SampleID","Location","Season", "Depth", "Dataset", #metadata
         "input","filtered", "merged", "nochim", "Seq.prop", #dada2 
         "Chloroplast","Mitochondria", "Archaea", #taxa
         "Observed","Chao1","Shannon","InvSimpson") #alpha div

write.table(summary_table2, "./tables/16S_15_20_overview_table_2.txt" , sep = "\t", quote = F)

####################################
#Plot rarefaction
####################################

#https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

library(iNEXT)
ps0.iNEXT <- iNEXT(as.data.frame(otu_table(ps0)), q=0, datatype="abundance", knots = 40)

ggiNEXT(ps0.iNEXT, type = 1)
ggsave("./output_graphs/16S_15_20_rarefactions.pdf", plot=last_plot())

ggiNEXT(ps0.iNEXT, type = 1, facet.var="site")

ggiNEXT(ps0.iNEXT, type = 2)
ggsave("./output_graphs/16S_15_20_rarefactions_coverage.pdf", plot=last_plot())


##################################
#Filtering
###################################

#Taxonomic filtering
####################

# Create table, number of features for each phyla
table_phyla <- table(tax_table(ps0)[, "Phylum"], exclude = NULL)
write.table(table_phyla, "./tables/TablePhyla.csv", quote = F)
write.table(table_phyla, "./tables/TablePhyla.txt", sep = "\t", quote = F)


#Remove features with a Phylum NA
phy_obj2 <- subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
phy_obj2


#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
# Compute prevalence of each feature, store as data.frame
prev = apply(X = otu_table(phy_obj2),
             MARGIN = ifelse(taxa_are_rows(phy_obj2), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(phy_obj2),
                    tax_table(phy_obj2))

prev_sum_phylum <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion=100*Abundance/sum(Abundance))
prev_sum_phylum
write.table(prev_sum_phylum, "./tables/prev_sum_phylum.txt", sep = "\t", quote = F)


#Define phyla to filter (observed only one feature)

fliterPhyla = c("Zixibacteria", "WS2", "LCP-89", "Modulibacteria", "Fermentibacterota", "Abditibacteriota", "Armatimonadota", "Halanaerobiaeota", "Acetothermia", "Caldisericota", "Deferrisomatota")
#filter phyla
phy_obj3 = subset_taxa (phy_obj2, !Phylum %in% fliterPhyla)
phy_obj3

#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
prevdf3 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj3, "Phylum"))
ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(phy_obj3),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") 

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts.pdf", last_plot())


#Prevalence filtering
########################

#Define prevalence treshold as 1 sample 
prevdf2 <- prevdf3 %>% filter (Prevalence>1)
prevdf2

#Remove phyla that have less than 0.01% of seq
# First agglomerate by phylum to define the phyla you want to keep.
physeqPhylum = tax_glom(phy_obj3, "Phylum")
physeqPhylumRA = transform_sample_counts(physeqPhylum, function(x) x/sum(x))
physeqPhylumRAF = filter_taxa(physeqPhylumRA, function(x) mean(x) > 0.0001, TRUE)

# Define the vector of phyla names that are still present after your filter
keepPhyla = get_taxa_unique(physeqPhylumRAF, "Phylum")
# Like in the first question, subset to just the phyla that you want, using the original phyloseq object
phy_obj4 = subset_taxa(phy_obj3, Phylum  %in% keepPhyla) 
phy_obj4
prevdf4 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj4, "Phylum"))
ggplot(prevdf4, aes(TotalAbundance, Prevalence / nsamples(phy_obj4),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#Remove phyla that have less than 0.1% of seq
# First agglomerate by phylum to define the phyla you want to keep.
physeqPhylum = tax_glom(phy_obj3, "Phylum")
physeqPhylumRA = transform_sample_counts(physeqPhylum, function(x) x/sum(x))
physeqPhylumRAF = filter_taxa(physeqPhylumRA, function(x) mean(x) > 0.001, TRUE)

# Define the vector of phyla names that are still present after your filter
keepPhyla = get_taxa_unique(physeqPhylumRAF, "Phylum")
# Like in the first question, subset to just the phyla that you want, using the original phyloseq object
phy_obj5 = subset_taxa(phy_obj3, Phylum  %in% keepPhyla) 
phy_obj5
prevdf5 = subset(prevdf, Phylum %in% get_taxa_unique(phy_obj5, "Phylum"))
ggplot(prevdf5, aes(TotalAbundance, Prevalence / nsamples(phy_obj5),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


#####################################
# Explore results on phylum, class and family level
#####################################

#Compare on Phylum level
prevdf_sum_phylum <- plyr::ddply(prevdf3, "Phylum", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
prevdf_sum_phylum4 <- plyr::ddply(prevdf4, "Phylum", function(prevdf4){cbind(Samples=mean(prevdf4$Prevalence),Abundance=sum(prevdf4$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
prevdf_sum_phylum5 <- plyr::ddply(prevdf5, "Phylum", function(prevdf5){cbind(Samples=mean(prevdf5$Prevalence),Abundance=sum(prevdf5$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#Explore results on class level
prevdf_sum_class <- plyr::ddply(prevdf3, "Class", function(prevdf3){cbind(Samples=mean(prevdf3$Prevalence),Abundance=sum(prevdf3$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
prevdf_sum_class4 <- plyr::ddply(prevdf4, "Class", function(prevdf4){cbind(Samples=mean(prevdf4$Prevalence),Abundance=sum(prevdf4$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))
prevdf_sum_class5 <- plyr::ddply(prevdf5, "Class", function(prevdf5){cbind(Samples=mean(prevdf5$Prevalence),Abundance=sum(prevdf5$TotalAbundance))})%>%
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
phy_obj4ra = transform_sample_counts(phy_obj4, function (otu) {otu/sum(otu)})
phy_obj5ra = transform_sample_counts(phy_obj5, function (otu) {otu/sum(otu)})

#Melt data
phy_obj3ra.melt <- psmelt(phy_obj3ra)
phy_obj4ra.melt <- psmelt(phy_obj4ra)
phy_obj5ra.melt <- psmelt(phy_obj5ra)
write.table(phy_obj5ra.melt, "./tables/phy_obj5ra.melt.txt")

#Calculate abundance for each taxa
phy_obj5_melt.agg.genus <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family+Genus, phy_obj5ra.melt,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
phy_obj5_melt.agg.genus$Abundance <- phy_obj5_melt.agg.genus$Abundance.sum*100
phy_obj5_melt.agg.genus<- phy_obj5_melt.agg.genus[phy_obj5_melt.agg.genus$Abundance.sum>0,]


##calculate abundance for each Family
phy_obj5_melt.agg.family <- as.data.frame(as.list(aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, phy_obj5ra.melt,
                                                   FUN = function(x) c(sum = sum(x), count=length(x)))))

#phy_obj5_melt.agg.family <- aggregate(Abundance~Location+Depth+Season+Date+Class+Order+Family, phy_obj5ra.melt, FUN="sum")
phy_obj5_melt.agg.family$Abundance <- phy_obj5_melt.agg.family$Abundance.sum*100
phy_obj5_melt.agg.family<- phy_obj5_melt.agg.family[phy_obj5_melt.agg.family$Abundance.sum>0,]
#aggregate different samplings in one season
#phy_obj5_melt.agg.family<- aggregate(Abundance~Location+Depth+Season+Class+Order+Family, phy_obj5_melt.agg.family, FUN="mean")


Campyl.melt <- subset(phy_obj5_melt.agg.family , subset = Class == "Campylobacteria")

Vibrio<-subset(phy_obj5_melt.agg.genus, subset = Genus == "Vibrio")
write.csv(Vibrio)

###########################
#Taxonomic compositions on class level
###########################
#Calculate abundance for each Class
phy_obj5_melt.agg.class <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj5ra.melt, FUN="sum")
phy_obj5_melt.agg.class$Abundance <- phy_obj5_melt.agg.class$Abundance*100
phy_obj5_melt.agg.class<- phy_obj5_melt.agg.class[phy_obj5_melt.agg.class$Abundance>0,]

#remove below 3% ra
threshold<- 3
phy_obj5_melt.agg.class$Class <- as.character(phy_obj5_melt.agg.class$Class)
taxa_classes <- unique(phy_obj5_melt.agg.class$Class[!phy_obj5_melt.agg.class$Abundance<threshold])
phy_obj5_melt.agg.class$Class[phy_obj5_melt.agg.class$Abundance<threshold] <- "Other taxa"
phy_obj5_melt.agg.class$Class <- factor(phy_obj5_melt.agg.class$Class,
                                        levels=c(taxa_classes,"Other taxa"))
phy_obj5_melt.agg.class.table <- aggregate (Abundance~Location+Season+Date+Depth+Class, phy_obj5_melt.agg.class, FUN="sum")
write.csv(phy_obj5_melt.agg.class.table, "./tables/Class.csv")


#Plot Class(Location)
phy_obj5_melt.agg.class$Location = factor(phy_obj5_melt.agg.class$Location, levels=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggplot(phy_obj5_melt.agg.class, aes(x = Abundance, y = Season, fill = Class))+
  facet_grid(Location~., space= "fixed")+
  geom_bar(stat = "identity", position="fill") +
  scale_y_discrete(limits = c('autumn', 'winter','summer','spring'))
ggsave("./output_graphs/Class5_location.pdf", last_plot())


#Plot Class(Season)
colourCount = length(unique(phy_obj5_melt.agg.class$Class))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
phy_obj5_melt.agg.class$Season_f = factor(phy_obj5_melt.agg.class$Season, levels=c('autumn','winter', 'summer','spring'))
ggplot(phy_obj5_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Class5_season_depth.pdf", last_plot())


#Plot Class(Date)
phy_obj5_melt.agg.class$Season = factor(phy_obj5_melt.agg.class$Date, levels=c( '15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
ggplot(phy_obj5_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Class5_Date.pdf", last_plot())


#Subsets by taxonomy - Campylobacteria (position stack/positiona fill)
Campyl.melt$Date_f = factor(Campyl.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
ggplot(Campyl.melt, aes(x = Abundance, y = Location, fill = Family))+
  facet_grid(Date_f~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Campylobacteria.pdf", last_plot())

phy_obj5_melt.agg.class$Season = factor(phy_obj5_melt.agg.class$Date, levels=c( '15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))
ggplot(phy_obj5_melt.agg.class, aes(x = Abundance, y = Location, fill = Class))+
  facet_grid(Season~Depth, space= "fixed")+
  geom_bar(stat = "identity", position="fill")+
  scale_fill_manual(values=getPalette(colourCount))+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Class5_Date.pdf", last_plot())


#Subsets by taxonomy - Alphaproteobacteria 
Alpha.melt <- subset(x=phy_obj5_melt.agg.family , subset = Class == "Alphaproteobacteria")
Alpha.melt$Date_f = factor(Alpha.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

colourCount = length(unique(Alpha.melt$Order))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

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

#Subsets by taxonomy - Gammaproteobacteria 
Gamma.melt <- subset(x=phy_obj5_melt.agg.family , subset = Class == "Gammaproteobacteria")
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


#Subsets by taxonomy - Bacteroidia
Bacter.melt <- subset(phy_obj5_melt.agg.family , subset = Class == "Bacteroidia")
Bacter.melt$Date_f = factor(Bacter.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Bacter.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

#Subsets by taxonomy - Bacilli
Bacilli.melt <- subset(phy_obj5_melt.agg.family , subset = Class == "Bacilli")
Bacilli.melt$Date_f = factor(Bacilli.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Bacilli.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Bacteroidia.pdf", last_plot())

#Subsets by taxonomy - Clostridia
Clost.melt <- subset(phy_obj5_melt.agg.family , subset = Class == "Clostridia")
Clost.melt$Date_f = factor(Clost.melt$Date, levels=c('15.10.2019','25.11.2015','4.12.2018','4.02.2016','24.03.2015','25.04.2019','18.06.2019','15.07.2015'))

ggplot(Clost.melt, aes(x = Abundance, y = Location, fill = Order))+
  facet_grid(Date_f~Depth, space= "fixed")+
  scale_fill_manual(values=getPalette(colourCount))+
  geom_bar(stat = "identity", position="stack")+
  scale_y_discrete(limits=c('SM-Outfall','OS-Marine','NS-Marine','R-Estuary-2','R-Estuary-1','R-Mouth'))
ggsave("./output_graphs/Clostridia.pdf", last_plot())


#Plot nMDS
phy_obj5ra.surface.2020 = subset_samples(phy_obj5ra, Depth == "surface" & Dataset=="2020")
phy_obj5.nmds <- ordinate(phy_obj5ra.surface.2020, method = "NMDS", 
                          distance = "jsd")

plot_ordination(phy_obj5ra.surface.2020, phy_obj5.nmds, color = "Season")+
  stat_ellipse(type="norm")+
  geom_point(aes(shape = Location),color ="black", size = 5)+
  geom_point(aes(shape = Location), size = 4)+
  coord_fixed()
ggsave("./output_graphs/nMDS_2020_surface.pdf", last_plot())

phy_obj5ra.bottom.2020 = subset_samples(phy_obj5ra, Depth == "bottom" & Dataset=="2020")
phy_obj5.nmds <- ordinate(phy_obj5ra.bottom.2020, method = "NMDS", 
                          distance = "jsd")
plot_ordination(phy_obj5ra.bottom.2020, phy_obj5.nmds, color = "Season")+
  stat_ellipse(type="norm")+
  geom_point(aes(shape = Location),color ="black", size = 5)+
  geom_point(aes(shape = Location), size = 4)+
  coord_fixed()

ggsave("./output_graphs/nMDS_bottom_2020.pdf", last_plot())



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

head(sample_data(phy_obj5)$Location)

phy_obj5.ds = subset_samples(phy_obj3, Location%in%c("R-Estuary-1", "NS-Marine")&Depth == "surface"&Dataset==2020)
head(sample_data(phy_obj5.ds)$Location)
head(sample_data(phy_obj5.ds)$Depth)
head(sample_data(phy_obj5.ds)$Dataset)

diagdds = phyloseq_to_deseq2(phy_obj5.ds, ~ Location)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#Investigate test results table
resultsNames(diagdds)

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy_obj5.ds)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#Plot summary of the results=ASVs that were significantly different between the two Locations


# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  scale_fill_manual(values=getPalette(colourCount))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Family order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  scale_fill_manual(values=getPalette(colourCount))+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


##Count how many NAs increase and decrease in abundance
sigtab.pos = sigtab[(sigtab$log2FoldChange >=0),]
pos.row.has.na <- apply(sigtab.pos, 1, function(x){any(is.na(x))})
sum(pos.row.has.na)

sigtab.neg = sigtab[(sigtab$log2FoldChange <=0),]
neg.row.has.na <- apply(sigtab.neg, 1, function(x){any(is.na(x))})
sum(neg.row.has.na)

ggsave("./output_graphs/DESeq2.pdf", last_plot())


#https://micca.readthedocs.io/en/latest/phyloseq.html
