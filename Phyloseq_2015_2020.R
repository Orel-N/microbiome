##############################
#Phyloseq Microbiome 2015_2020
##############################
#instal phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq","dplyr", "rstatix", "iNEXT"))

library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(RColorBrewer)
library(iNEXT)
library(DESeq2)
library(rstatix)
library(vegan)

theme_set(theme_bw())

##################################
#Import dada2 output into phyloseq
##################################

#Import Sample Data - "metadata"
meta <- read.csv("./data/Metadata_2015_2020.csv", h=T, sep = ",")

#Import "ASV table" (matrix)
seqtab.nochim <- read.csv("./dada2/dada2_seqtab_nochim2_2.txt", h=T, sep="\t")

#Import taxonomy table (matrix)
taxa <- as.matrix(read.csv("./dada2/dada2_taxonomy_table_2.txt", h=T,sep = "\t"))

#Edit metadata
# Checking data structure
summary(meta)
str(meta)

rownames(meta)<- paste("X",meta$SampleID, sep ="")
meta
meta[,'DIN']=round(meta[,'DIN'],1)
meta[,'Depth_m']=round(meta[,'Depth_m'],0)
meta[,'Coliform_bacteria']=round(meta[,'Coliform_bacteria'],0)
meta[,'Cyanobacteria_STDEV']=round(meta[,'Cyanobacteria_STDEV'],0)

meta$Dataset <- ifelse(meta$SampleID < 12, "2015", "2020")
                                  
#Check order
all.equal(rownames(seqtab.nochim), rownames(taxa))

#Create a phyloseq object from the OTU table/ASV table and taxonomy assigned by DADA2
ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), 
              sample_data(meta), tax_table(taxa))
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
pruned_seq_sums$Row.names <- paste(row.names(pruned_seq_sums))
pruned_seq_sums

#remove unclassified on phylum level, chloroplast, Mitochondrial and Archaeal sequence variants
ps0 <- subset_taxa(ps, !Kingdom %in% c("Eukaryota") & !is.na(Phylum) & !Phylum %in% c("", "uncharacterized") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") & !Kingdom %in% c("Archaea"))
ps0

#Create separate phyloseq object with seq assigned to Chloroplast
ps0_chl <- subset_taxa(ps, Order == "Chloroplast")
ps0_chl


###########################################
#Alpha diversity
###########################################

#Dataset 2015 and 2020 together
####################################

#Create Alphadiversity table
ps0_alpha  <- estimate_richness(ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
alpha_df_15.20 <- data.frame(sample_data(ps0), ps0_alpha)
write.table(alpha_df_15.20, "./tables/AlphaDiversity15_20.txt")
ps0_alpha$Row.names<-rownames(ps0_alpha)

#Edit table
alpha_df_15.20$Season <- factor(alpha_df_15.20$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_15.20$Location <- factor(alpha_df_15.20$Location, levels = c("R-Mouth", "R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_15.20$Depth <- factor(alpha_df_15.20$Depth, levels = c("bottom", "surface"))
alpha_df_15.20$Pollution <- ifelse(alpha_df_15.20$Location == "OS-Marine", "Unpolluted", 
                                   ifelse(alpha_df_15.20$Location == "NS-Marine", "Unpolluted", "Polluted"))
#Plot Chao1 - separate for depth
Chao.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.p
ggsave("./output_graphs/Chao_15_20.pdf")

#Plot Shannon d.i. - separate for depth
Shannon.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.p
ggsave("./output_graphs/Shannon_15_20.pdf")

#Plot InvSimpson - separate for depth
InvSimpson.p <- ggplot(data = alpha_df_15.20) +
  geom_boxplot(mapping=aes(x=Location, y=InvSimpson, color=Location)) +
  geom_point(mapping=aes(x=Location, y=InvSimpson, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
InvSimpson.p
ggsave("./output_graphs/InvSimpsin_15_20.pdf")

ggarrange(Chao.p, Shannon.p, InvSimpson.p, nrow=1, common.legend = TRUE)

###Statistic

#Shapiro-Wilk's test
hist(alpha_df_15.20$Chao1)
qqnorm(alpha_df_15.20$Chao1)
qqline(alpha_df_15.20$Chao1)
hist(alpha_df_15.20$Shannon)
qqnorm(alpha_df_15.20$Shannon)
qqline(alpha_df_15.20$Shannon)

alpha_df_15.20 %>% shapiro_test(Chao1, Shannon, InvSimpson)

anova.Ch.15.20 <- aov(alpha_df_15.20$Chao1 ~ Dataset+Depth+Season+Pollution+Location, data=alpha_df_15.20)
summary(anova.Ch.15.20)
TukeyHSD(anova.Ch.15.20)

kruskal.test(alpha_df_15.20$Shannon ~ Season, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$Shannon ~ Depth, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$Shannon ~ Location, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$Shannon ~ Dataset, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$Shannon ~ Pollution, data = alpha_df_15.20)

kruskal.test(alpha_df_15.20$InvSimpson ~ Season, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$InvSimpson ~ Depth, data = alpha_df_15.20) 
kruskal.test(alpha_df_15.20$InvSimpson ~ Location, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$InvSimpson ~ Dataset, data = alpha_df_15.20)
kruskal.test(alpha_df_15.20$InvSimpson ~ Pollution, data = alpha_df_15.20)

pairwise.wilcox.test(alpha_df_15.20$InvSimpson, alpha_df_15.20$Season,
                     p.adjust.method = "BH")

#Comparing just surface samples from both datasets
alpha_df_15.20.s <- filter(alpha_df_15.20, Depth == "surface")
alpha_df_15.20.s %>% shapiro_test(Chao1, Shannon, InvSimpson)
anova.Ch.15.20.s <- aov(alpha_df_15.20.s$Chao1 ~ Dataset+Season+Pollution+Location, data=alpha_df_15.20.s)
summary(anova.Ch.15.20.s)
anova.Sh.15.20.s <- aov(alpha_df_15.20.s$Shannon ~ Dataset+Season+Pollution+Location, data=alpha_df_15.20.s)
summary(anova.Sh.15.20.s)
kruskal.test(alpha_df_15.20.s$InvSimpson ~ Dataset, data = alpha_df_15.20.s)

#Dataset 2020
####################################

ps0_20 <- subset_samples(ps0, Dataset == "2020")

#Create Alphadiversity table 2020
alpha_df_2 <- estimate_richness(ps0_20, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
alpha_df_20 <- data.frame(sample_data(ps0_20), alpha_df_2)
write.table(alpha_df_20, "./tables/AlphaDiversity_20.txt")


#Edit table 2020
alpha_df_20$Season <- factor(alpha_df_20$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_20$Location <- factor(alpha_df_20$Location, levels = c("R-Mouth", "R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_20$Depth <- factor(alpha_df_20$Depth, levels = c("bottom", "surface"))
alpha_df_20$Pollution <- ifelse(alpha_df_20$Location == "OS-Marine", "Unpolluted", 
                                ifelse(alpha_df_20$Location == "NS-Marine", "Unpolluted", "Polluted"))

#Plot Chao1 - separate for depth 2020
Chao.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.20.p
ggsave("./output_graphs/Chao_20.pdf")

#Plot Shannon d.i. - separate for depth 2020
Shannon.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.20.p
ggsave("./output_graphs/Shannon_20.pdf")

#Plot InvSimpson - separate for depth 2020
InvSimpson.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=InvSimpson, color=Location)) +
  geom_point(mapping=aes(x=Location, y=InvSimpson, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
InvSimpson.20.p
ggsave("./output_graphs/InvSimpsin_20.pdf")

alpha.20 <- ggarrange(Chao.20.p, Shannon.20.p, InvSimpson.20.p, nrow=1, common.legend = TRUE)
alpha.20

###Statistic 2020
#Subsetting data 2020 - bottom vs.surface
alpha_df_20b <- filter(alpha_df_20, Depth == "bottom")
alpha_df_20s <- filter(alpha_df_20, Depth == "surface")

#Shapiro-Wilk's test
hist(alpha_df_20$Chao1)
qqnorm(alpha_df_20$Chao1)
qqline(alpha_df_20$Chao1)
hist(alpha_df_20$Shannon)
qqnorm(alpha_df_20$Shannon)
qqline(alpha_df_20$Shannon)

alpha_df_20 %>% shapiro_test(Chao1, Shannon, InvSimpson)
alpha_df_20b %>% shapiro_test(Chao1, Shannon, InvSimpson)
alpha_df_20s %>% shapiro_test(Chao1, Shannon, InvSimpson)

#ANOVA or Kruskal-Wallis test
anova.Ch.20 <- aov(alpha_df_20$Chao1 ~ Depth+Season+Pollution+Location, data=alpha_df_20)
summary(anova.Ch.20)
TukeyHSD(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20b$Chao1 ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20s$Chao1 ~ Season+Pollution+Location, data=alpha_df_20s)
summary(anova.Ch.20)

kruskal.test(alpha_df_20$Shannon ~ Season, data = alpha_df_20) 
kruskal.test(alpha_df_20$Shannon ~ Depth, data = alpha_df_20) 
kruskal.test(alpha_df_20$Shannon ~ Location, data = alpha_df_20)
kruskal.test(alpha_df_20$Shannon ~ Pollution, data = alpha_df_20)
anova.Ch.20 <- aov(alpha_df_20b$Shannon ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.Ch.20)
anova.Ch.20 <- aov(alpha_df_20s$Shannon ~ Season+Pollution+Location, data=alpha_df_20s)
summary(anova.Ch.20)
TukeyHSD(anova.Ch.20)

anova.InvS.20 <- aov(alpha_df_20$InvSimpson ~ Depth+Season+Pollution+Location, data=alpha_df_20)
summary(anova.InvS.20)
anova.InvS.20b <- aov(alpha_df_20b$InvSimpson ~ Season+Pollution+Location, data=alpha_df_20b)
summary(anova.InvS.20b)
anova.InvS.20s <- aov(alpha_df_20s$InvSimpson ~ Season+Pollution+Location, data=alpha_df_20s)
summary(anova.InvS.20s)



#Dataset 2015
####################################

ps0_15 <- subset_samples(ps0, Dataset == "2015")


#Create Alphadiversity table 2015
alpha_df_3 <- estimate_richness(ps0_15, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
alpha_df_15 <- data.frame(sample_data(ps0_15), alpha_df_3)
write.table(alpha_df_15, "./tables/AlphaDiversity_15.txt")


#Edit table 2015
alpha_df_15$Season <- factor(alpha_df_15$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_15$Location <- factor(alpha_df_15$Location, levels = c("R-Mouth", "R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_15$Depth <- factor(alpha_df_15$Depth, levels = c("bottom", "surface"))
alpha_df_15$Pollution <- ifelse(alpha_df_15$Location == "OS-Marine", "Unpolluted", 
                                ifelse(alpha_df_15$Location == "NS-Marine", "Unpolluted", "Polluted"))

#Plot Chao1 - separate for depth 2015
Chao.15.p <- ggplot(data = alpha_df_15) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.15.p
ggsave("./output_graphs/Chao_15.pdf")

#Plot Shannon d.i. - separate for depth 2015
Shannon.15.p <- ggplot(data = alpha_df_15) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.15.p
ggsave("./output_graphs/Shannon_15.pdf")

#Plot InvSimpson - separate for depth 2015
InvSimpson.15.p <- ggplot(data = alpha_df_15) +
  geom_boxplot(mapping=aes(x=Location, y=InvSimpson, color=Location)) +
  geom_point(mapping=aes(x=Location, y=InvSimpson, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
InvSimpson.15.p
ggsave("./output_graphs/InvSimpsin_15.pdf")

alpha.15 <- ggarrange(Chao.15.p, Shannon.15.p, InvSimpson.15.p, nrow=1, common.legend = TRUE)
alpha.15

###Statistic 2015

#Shapiro-Wilk's test
hist(alpha_df_15$Chao1)
qqnorm(alpha_df_15$Chao1)
qqline(alpha_df_15$Chao1)
hist(alpha_df_15$Shannon)
qqnorm(alpha_df_15$Shannon)
qqline(alpha_df_15$Shannon)

alpha_df_15 %>% shapiro_test(Chao1, Shannon, InvSimpson)

#ANOVA or Kruskal-Wallis test
anova.Ch.15 <- aov(alpha_df_15$Chao1 ~ Pollution+Location+Season, data=alpha_df_15)
summary(anova.Ch.15)
TukeyHSD(anova.Ch.15)

anova.Sh.15 <- aov(alpha_df_15$Shannon ~ Pollution+Location+Season, data=alpha_df_15)
summary(anova.Sh.15)
TukeyHSD(anova.Sh.15)

kruskal.test(alpha_df_15$InvSimpson ~ Season, data = alpha_df_15) 
kruskal.test(alpha_df_15$InvSimpson ~ Location, data = alpha_df_15)
kruskal.test(alpha_df_15$InvSimpson ~ Pollution, data = alpha_df_15)


##################
#Overview table

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
ps0.iNEXT <- iNEXT(as.data.frame(otu_table(ps0)), q=0, datatype="abundance", knots = 40)

#Plot all together
ggiNEXT(ps0.iNEXT, type = 1)
ggsave("./output_graphs/16S_15_20_rarefactions.pdf", plot=last_plot())

ggiNEXT(ps0.iNEXT, type = 1, facet.var="site")
ggiNEXT(ps0.iNEXT, type = 2)
ggsave("./output_graphs/16S_15_20_rarefactions_coverage.pdf", plot=last_plot())


#Plot seperatly Depth~Location
ps0.iNEXT.rare <-fortify(ps0.iNEXT, type=1)
ps0.iNEXT.meta <- as(sample_data(ps0), "data.frame")
ps0.iNEXT.meta$site <- rownames(ps0.iNEXT.meta)

ps0.iNEXT.rare$Location <- ps0.iNEXT.meta$Location[match(ps0.iNEXT.rare$site, ps0.iNEXT.meta$site)]
ps0.iNEXT.rare$Depth <- ps0.iNEXT.meta$Depth[match(ps0.iNEXT.rare$site, ps0.iNEXT.meta$site)]
ps0.iNEXT.rare$Season <- ps0.iNEXT.meta$Season[match(ps0.iNEXT.rare$site, ps0.iNEXT.meta$site)] 

ps0.iNEXT.rare.point <- ps0.iNEXT.rare[which(ps0.iNEXT.rare$method == "observed"),]
ps0.iNEXT.rare.line <- ps0.iNEXT.rare[which(ps0.iNEXT.rare$method != "observed"),]
ps0.iNEXT.rare.line$method <- factor (ps0.iNEXT.rare.line$method,
                                      c("interpolated", "extrapolated"),
                                      c("interpolation", "extrapolation"))
iNEXT.rare.line <- data.frame()
iNEXT.rare.point<- data.frame()

iNEXT.rare.line <- rbind(iNEXT.rare.line,ps0.iNEXT.rare.line)
iNEXT.rare.point <- rbind(iNEXT.rare.point,ps0.iNEXT.rare.point)

iNEXT.rare.line$Season <- factor(iNEXT.rare.line$Season, levels = c("autumn", "winter", "spring", "summer"))
iNEXT.rare.point$Season <- factor(iNEXT.rare.point$Season, levels = c("autumn", "winter", "spring", "summer"))

rare.p <- ggplot(iNEXT.rare.line, aes(x=x, y=y, shape = site))+
  geom_line(aes(linetype = method, colour = Season), lwd = 1, data= iNEXT.rare.line)+
  geom_point(shape = 21, size =3, colour = "black", data= iNEXT.rare.point)+
  labs(x = "Sample size", y = "Species richness")+
  xlim(0,3e5)+
  theme_classic(base_size = 12)+theme(legend.position="top")+
  facet_grid(Depth~Location, scales = "fixed",as.table = TRUE)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) + 
  coord_cartesian(clip="off")
rare.p

ggsave("./output_graphs/16S_15_20_rarefactions_location_depth.pdf", plot=last_plot())


##################################
#Filtering
###################################

#Taxonomic filtering
####################

# Create table, number of features for each phyla
table_phyla <- table(tax_table(ps0)[, "Phylum"], exclude = NULL)
write.table(table_phyla, "./tables/TablePhyla2.csv", quote = F)
write.table(table_phyla, "./tables/TablePhyla2.txt", sep = "\t", quote = F)

phy_obj2 <- ps0

#Explore feature prevalence in the dataset - number of samples in which a taxa appears at least once
# Compute prevalence of each feature, store as data.frame
prev = apply(X = otu_table(phy_obj2),
             MARGIN = ifelse(taxa_are_rows(phy_obj2), yes = 1, no = 2),
             FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prev,
                    TotalAbundance = taxa_sums(phy_obj2),
                    tax_table(phy_obj2))
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(phy_obj2),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts_BeforeF.pdf", last_plot())

prev_sum_phylum <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion=100*Abundance/sum(Abundance))
prev_sum_phylum
write.table(prev_sum_phylum, "./tables/prev_sum_phylum.txt", sep = "\t", quote = F)

#Define phyla to filter (observed only one feature)
fliterPhyla = c("Sumerlaeota", "Cloacimonadota", "Schekmanbacteria", "Hydrogenedentes", "Deinococcota", "Thermotogota","Calditrichota","Zixibacteria", "WS2", "LCP-89", "Modulibacteria", "Fermentibacterota", "Abditibacteriota", "Armatimonadota", "Halanaerobiaeota", "Acetothermia", "Caldisericota", "Deferrisomatota")
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

ggsave("./output_graphs/TaxaPrevalenceVsTotalCounts_AfterF.pdf", last_plot())

#Save filtered phyloseq object
saveRDS(phy_obj3, "./phyloseqFiltered.RDS")
saveRDS(prevdf3, "./prevdfFiltered.RDS")
saveRDS(ps0_chl, "./chloroplasts.RDS")



