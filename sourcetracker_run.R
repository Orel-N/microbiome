library("phyloseq")
library("dplyr")
library("reshape2")
library("tidyr")
library("ggplot2")

#source('./scripts/SourceTracker_EF.R')
source('./scripts/Color_palettes.R')
source('./scripts/SourceTracker.r')

phy_obj3<- readRDS("./data/phyloseqPrevFiltered.RDS")

phy_obj3<- subset_samples(phy_obj3, Dataset == "2020")
phy_obj3<- prune_taxa(taxa_sums(phy_obj3)>0,phy_obj3)

#ASV table
otus <- t(as(otu_table(phy_obj3, taxa_are_rows = FALSE),"matrix"))

#metadata
metadata <- as(sample_data(phy_obj3),"data.frame") 

# extract the source environments and source/sink indices
train.ix <- which(metadata$Category=='source')
test.ix <- which(metadata$Category=='sink')
# extract the source environments and source/sink indices
train.ix <- which(metadata$Category=='source')
test.ix <- which(metadata$Category=='sink')
envs <- metadata$Location

#####################################
#Run SourceTracker
#####################################
#no need to run this section, the results already in the data directory:
load("./data/ST_2020.RData")

# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(otus[train.ix,], envs[train.ix])
#alpha1 <- tune.results$best.alpha1
#alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st_no_season <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth=1000)

# Estimate source proportions in test data
ST_results_no_season <- predict.sourcetracker(st_no_season,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

# Estimate leave-one-out source proportions in training data 
ST_results.train_no_season <- predict.sourcetracker(st_no_season, alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

#####################################
#Plot results
#####################################

#make labels
labels <- sprintf('%s %s %s', metadata$Location, metadata$Season, metadata$Depth)

#rename the the results tables
rownames(ST_results_no_season$proportions) <- labels[test.ix]
rownames(ST_results.train_no_season$proportions) <- labels[train.ix]

#melt and combined the train and test datasets
results_melted <- melt(ST_results_no_season$proportions)
results.train_melted <- melt(ST_results.train_no_season$proportions)

#merge tables and define variables order
resluts_merged <- rbind(results_melted,results.train_melted) %>% 
  separate(Var1, c("Location", "Season", "Depth"), " ") %>% 
  mutate(Season = factor(Season, levels = c("winter","spring", "summer","autumn")),
         Depth = factor(Depth, levels = c("surface", "bottom")),
         Location = factor(Location, levels = c("R-Estuary-1","R-Estuary-2", "NS-Marine", "OS-Marine", "SM-Outfall")))

#plot
resluts_merged_plot <- ggplot() + 
  geom_bar(aes(y = value, x = Location, fill = Var2), colour = "black", 
           data = resluts_merged, stat="identity")+
  facet_grid(Depth~Season)+
  labs(y="Average source contribution")+
  scale_fill_manual(values = sample(tol21rainbow))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle=90),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

ggsave("./Figures/mic_sourcetracker.png", 
       plot = resluts_merged_plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Explore affiliation of specific taxa to the different sources
#####################################
#To get a matrix of ASV x Probability(ASV came from source) across all samples
ASV.source.prob <- t(apply(ST_results.train_no_season$full.results,c(2,3),mean))
ASV.source.prob <- as.data.frame(sweep(ASV.source.prob,1,rowSums(ASV.source.prob),'/'))

#Add correct rownames and columnnames and add taxaID
taxatable <- as(tax_table(phy_obj3),"matrix")
taxatable <- as.data.frame(taxatable)

row.names(ASV.source.prob) <- row.names(taxatable)
colnames(ASV.source.prob) <- c("WW","SM","Unknown")

taxa_sourcetracker <- cbind(taxatable,ASV.source.prob) #%>% melt() %>% filter(value >0)

#check which and how many ASVs were affiliated to each source
taxa_sourcetracker_WW <- taxa_sourcetracker %>% 
  select(Phylum,Class,Order,Family,Genus,Species,WW) %>% 
  filter(WW> 0.99) %>% 
  group_by(Phylum,Class,Order,Family) %>% 
  summarise(ASVs=n())

taxa_sourcetracker_SM <- taxa_sourcetracker %>% 
  select(Phylum,Class,Order,Family,Genus,Species,SM) %>% 
  filter(SM> 0.99) %>% 
  group_by(Phylum,Class,Order,Family) %>% 
  summarise(ASVs=n())

taxa_sourcetracker_ASVs<- merge(taxa_sourcetracker_WW,taxa_sourcetracker_SM, by =c("Phylum","Class","Order","Family"), all = TRUE) %>% 
                            plyr::rename(c("ASVs.x" = "WW","ASVs.y" = "SM"))

write.csv(taxa_sourcetracker_ASVs,"./Tables/taxa_sourcetracker_ASVs.csv")

#Select microbial indicators
Mic_Ind <- read.table("./data/Microbial_Indicators.txt", h=T, sep="\t")
sourcetracker_mic_indicators <- taxa_sourcetracker %>% filter(Family %in% c(Mic_Ind$Family))

taxa_sourcetracker_ASVs_indicators<- taxa_sourcetracker_ASVs %>% filter(Family %in% c(Mic_Ind$Family))

write.csv(taxa_sourcetracker_ASVs_indicators,"./Tables/taxa_sourcetracker_mic_ind.csv")
