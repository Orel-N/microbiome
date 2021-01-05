source('./SourceTracker.R')

library("phyloseq")

phy_obj3<- readRDS("./phyloseqFiltered.RDS")

#ASV table
otus <- as.data.frame(otu_table(phy_obj3, taxa_are_rows = TRUE))
otus <- t(as.matrix(otus))

#metadata
metadata <- sample_data(phy_obj3)

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$Category=='source')
test.ix <- which(metadata$Category=='sink')
#envs <- metadata$Location
envs<- gsub("-.*","",metadata$Location)
if(is.element('Location',colnames(metadata))) desc <- metadata$Location
if(is.element('Location',colnames(metadata))) depth <- metadata$Depth

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
st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth=1000)

save.image("./ST_results.RData")

# Estimate source proportions in test data
ST_results <- predict.sourcetracker(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

save.image("./ST_results.RData")

# Estimate leave-one-out source proportions in training data 
ST_results.train <- predict.sourcetracker(st, alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, full.results=TRUE)

save.image("./ST_results.RData")
