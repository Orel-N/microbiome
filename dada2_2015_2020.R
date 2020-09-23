#################################
#dada2 16S Microbiome 2015, 2020
#################################
#Run R
R

#load libraries and set random seed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
set.seed(123)

#Define the path for the fastq files
path <- "/DATA/mikrobiom/16S_2015_2020_DADA2/Clipped" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

filt_path <- file.path("Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)


#Inspect read quality profiles
plotQualityProfile(fnFs[1:51])
ggsave("qualplotFs_2015_2020.pdf", last_plot())
plotQualityProfile(fnRs[1:51])
ggsave("qualplotRs_2015_2020.pdf", last_plot())

# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#separate the different runs
#2015
fnFs_15<- sort(file.path("Clipped",paste(c(1:11), "_clip_R1.fastq", sep = "")))
fnRs_15<- sort(file.path("Clipped",paste(c(1:11), "_clip_R2.fastq", sep = ""))) 
filtFs_15<- sort(file.path("Filtered",paste(c(1:11), "_F_filt.fastq.gz", sep = "")))
filtRs_15<- sort(file.path("Filtered",paste(c(1:11), "_R_filt.fastq.gz", sep = ""))) 

fnFs_15_2<- sort(file.path("Clipped",paste(c(1:11), "_clip_R1.fastq", sep = "")))
fnRs_15_2<- sort(file.path("Clipped",paste(c(1:11), "_clip_R2.fastq", sep = ""))) 
filtFs_15_2<- sort(file.path("Filtered",paste(c(1:11), "_F_filt.fastq.gz", sep = "")))
filtRs_15_2<- sort(file.path("Filtered",paste(c(1:11), "_R_filt.fastq.gz", sep = "")))

#Filter and trim
out_15 <- filterAndTrim(fnFs_15, filtFs_15, fnRs_15, filtRs_15, truncLen=c(220,220),
                        maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)

out_15_2 <- filterAndTrim(fnFs_15_2, filtFs_15_2, fnRs_15_2, filtRs_15_2,truncLen=c(210,210),
                      maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)
# Learn errors 
errF_15 <- learnErrors(filtFs_15, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_15 <- learnErrors(filtRs_15, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)

errF_15_2 <- learnErrors(filtFs_15_2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_15_2 <- learnErrors(filtRs_15_2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)



# Sample Inference 
dadaFs_15 <- dada(filtFs_15, err=errF_15, multithread=TRUE, verbose = TRUE)
dadaRs_15 <- dada(filtRs_15, err=errR_15, multithread=TRUE, verbose = TRUE)

dadaFs_15_2 <- dada(filtFs_15_2, err=errF_15_2, multithread=TRUE, verbose = TRUE)
dadaRs_15_2 <- dada(filtRs_15_2, err=errR_15_2, multithread=TRUE, verbose = TRUE)

#Merge paired reads
mergers_15 <- mergePairs(dadaFs_15, filtFs_15, dadaRs_15, filtRs_15, verbose=TRUE, minOverlap = 10)

mergers_15_2 <- mergePairs(dadaFs_15_2, filtFs_15_2, dadaRs_15_2, filtRs_15_2, verbose=TRUE, minOverlap = 10)

#make sequence table and save it
seqtab_15 <- makeSequenceTable(mergers_15)
saveRDS(seqtab_15, "seqtab_15.rds")


#Plot error profiles
plot.errF_15 <- plotErrors(errF_15, nominalQ=TRUE)+ggtitle("Error rates F 2015")
ggsave("plotErrorsFs2015.pdf", plot.errF_15, device="pdf")
plot.errR_15 <- plotErrors(errR_15, nominalQ=TRUE)+ggtitle("Error rates R 2015")
ggsave("plotErrorsRs2015.pdf", plot.errR_15, device="pdf")

#Summary table
getN <- function(x) sum(getUniques(x))
track <- cbind(out_15, sapply(dadaFs_15, getN), sapply(dadaRs_15, getN), sapply(mergers_15, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)

track_2 <- cbind(out_15, sapply(dadaFs_15_2, getN), sapply(dadaRs_15_2, getN), sapply(mergers_15_2, getN))
colnames(track_2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track_2) <- sample.names
head(track_2)

#2020
fnFs_20<- sort(file.path("Clipped",paste(c(12:51), "_clip_R1.fastq", sep = "")))
fnRs_20<- sort(file.path("Clipped",paste(c(12:51), "_clip_R2.fastq", sep = ""))) 
filtFs_20<- sort(file.path("Filtered",paste(c(12:51), "_F_filt.fastq.gz", sep = "")))
filtRs_20<- sort(file.path("Filtered",paste(c(12:51), "_R_filt.fastq.gz", sep = ""))) 

fnFs_20_2<- sort(file.path("Clipped",paste(c(12:51), "_clip_R1.fastq", sep = "")))
fnRs_20_2<- sort(file.path("Clipped",paste(c(12:51), "_clip_R2.fastq", sep = ""))) 
filtFs_20_2<- sort(file.path("Filtered",paste(c(12:51), "_F_filt.fastq.gz", sep = "")))
filtRs_20_2<- sort(file.path("Filtered",paste(c(12:51), "_R_filt.fastq.gz", sep = "")))

#Filter and trim
out_20 <- filterAndTrim(fnFs_20, filtFs_20, fnRs_20, filtRs_20, truncLen=c(240,240),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)

out_20_2 <- filterAndTrim(fnFs_20_2, filtFs_20_2, fnRs_20_2, filtRs_20_2,truncLen=c(230,220),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE)
# Learn errors 
errF_20 <- learnErrors(filtFs_20, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_20 <- learnErrors(filtRs_20, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)

errF_20_2 <- learnErrors(filtFs_20_2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_20_2 <- learnErrors(filtRs_20_2, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)



# Sample Inference 
dadaFs_20 <- dada(filtFs_20, err=errF_20, multithread=TRUE, verbose = TRUE)
dadaRs_20 <- dada(filtRs_20, err=errR_20, multithread=TRUE, verbose = TRUE)

dadaFs_20_2 <- dada(filtFs_20_2, err=errF_20_2, multithread=TRUE, verbose = TRUE)
dadaRs_20_2 <- dada(filtRs_20_2, err=errR_20_2, multithread=TRUE, verbose = TRUE)

#Merge paired reads
mergers_20 <- mergePairs(dadaFs_20, filtFs_20, dadaRs_20, filtRs_20, verbose=TRUE, minOverlap = 10)

mergers_20_2 <- mergePairs(dadaFs_20_2, filtFs_20_2, dadaRs_20_2, filtRs_20_2, verbose=TRUE, minOverlap = 10)

#make sequence table and save it
seqtab_20_2 <- makeSequenceTable(mergers_20_2)
saveRDS(seqtab_20_2, "seqtab_20_2.rds")

#Plot error profiles
plot.errF_15 <- plotErrors(errF_15, nominalQ=TRUE)+ggtitle("Error rates F 2015")
ggsave("plotErrorsFs2015.pdf", plot.errF_15, device="pdf")
plot.errR_15 <- plotErrors(errR_15, nominalQ=TRUE)+ggtitle("Error rates R 2015")
ggsave("plotErrorsRs2015.pdf", plot.errR_15, device="pdf")
plot.errF_15_2 <- plotErrors(errF_15_2, nominalQ=TRUE)+ggtitle("Error rates F 2015_2")
ggsave("plotErrorsFs2015_2.pdf", plot.errF_15_2, device="pdf")
plot.errR_15_2 <- plotErrors(errR_15_2, nominalQ=TRUE)+ggtitle("Error rates R 2015_2")
ggsave("plotErrorsRs2015_2.pdf", plot.errR_15_2, device="pdf")

plot.errF_20 <- plotErrors(errF_20, nominalQ=TRUE)+ggtitle("Error rates F 2020")
ggsave("plotErrorsFs2020.pdf", plot.errF_20, device="pdf")
plot.errR_20 <- plotErrors(errR_20, nominalQ=TRUE)+ggtitle("Error rates R 2020")
ggsave("plotErrorsRs2020.pdf", plot.errR_20, device="pdf")
plot.errF_20_2 <- plotErrors(errF_20_2, nominalQ=TRUE)+ggtitle("Error rates F 2020_2")
ggsave("plotErrorsFs2020_2.pdf", plot.errF_20_2, device="pdf")
plot.errR_20_2 <- plotErrors(errR_20_2, nominalQ=TRUE)+ggtitle("Error rates R 2020_2")
ggsave("plotErrorsRs2020_2.pdf", plot.errR_20_2, device="pdf")

#Summary table
sample.names_15 <- sapply(strsplit(basename(fnFs_15), "_"), `[`, 1)
sample.names_15

sample.names_20 <- sapply(strsplit(basename(fnFs_20), "_"), `[`, 1)
sample.names_20

sample.names_15<-paste(1:11)
getN <- function(x) sum(getUniques(x))
track_15 <- cbind(out_15, sapply(dadaFs_15, getN), sapply(dadaRs_15, getN), sapply(mergers_15, getN))
colnames(track_15) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track_15) <- sample.names_15
head(track_15)


track_15_2 <- cbind(out_15_2, sapply(dadaFs_15_2, getN), sapply(dadaRs_15_2, getN), sapply(mergers_15_2, getN))
colnames(track_15_2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track_15_2) <- sample.names_15
head(track_15_2)

track_20 <- cbind(out_20, sapply(dadaFs_20, getN), sapply(dadaRs_20, getN), sapply(mergers_20, getN))
colnames(track_20) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track_20) <- sample.names_20
head(track_20)

track_20_2 <- cbind(out_20_2, sapply(dadaFs_20_2, getN), sapply(dadaRs_20_2, getN), sapply(mergers_20_2, getN))
colnames(track_20_2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track_20_2) <- sample.names_20
head(track_20_2)

write.table(track_15, "Overview_15.txt" , sep = "\t", quote = F)
write.table(track_20, "Overview_20.txt" , sep = "\t", quote = F)
write.table(track_15_2, "Overview_15_2.txt" , sep = "\t", quote = F)
write.table(track_20_2, "Overview_20_2.txt" , sep = "\t", quote = F)


#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1=makeSequenceTable(mergers_15),
                             table2=makeSequenceTable(mergers_20_2))


#Combine together sequences that are identical (run as separate script) 
##run as a seperate script
load("16S_2015_2020_DADA2.Rdata")
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)
dim(seqtab1)

save.image("16S_2015_2020_DADA2_2.Rdata")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

seqtab.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab1)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V34 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(380:430) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))
rownames(seqtab.nochim2) <- sample.names
saveRDS(seqtab.nochim2, file.path("./dada2", "seqtab15_20.rds"))

#Generate overview table - all together
total <- rbind(track_15, track_20_2)
total <- as.data.frame(total)


total$nochim<-colSums(seqtab.nochim)

write.table(total, "./data/Overview_dada2_15_20.txt" , sep = "\t", quote = F)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../16S_2015_DADA2/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "../16S_2015_DADA2/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)
saveRDS(taxa, file.path("./data", "taxa15_20.rds"))

##Write output
write.table(taxa, "dada2/dada2_taxonomy_table.txt", sep = "\t", quote = F)
write.table(t(seqtab.nochim2), "dada2/dada2_seqtab_nochim2.txt", quote = F, sep = "\t")

save.image("16S_2015_2020_DADA2_2.Rdata")

#quit R
q()