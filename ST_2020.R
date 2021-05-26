source('./SourceTracker.R')

library("phyloseq")


phy_obj3<- readRDS("./phyloseqPrevFiltered.RDS")
phy_obj3<- subset_samples(phy_obj3, Dataset == "2020")
phy_obj3<- prune_taxa(taxa_sums(phy_obj3)>0,phy_obj3)

#ASV table
otus <- t(as(otu_table(phy_obj3, taxa_are_rows = FALSE),"matrix"))

#metadata
metadata <- as(sample_data(phy_obj3),"data.frame") 

# extract the source environments and source/sink indices
train.ix <- which(metadata$Category=='source')
test.ix <- which(metadata$Category=='sink')
envs <- paste(metadata$Location,metadata$Season, sep ="_")

#####################################
#Run SourceTracker
#####################################
# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(otus[train.ix,], envs[train.ix])
#alpha1 <- tune.results$best.alpha1
#alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth=1000)

# Estimate source proportions in test data
ST_results <- predict.sourcetracker(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

# Estimate leave-one-out source proportions in training data 
ST_results.train <- predict.sourcetracker(st, alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

save.image("ST_2020.RData")
