#########################
#dada2 16S V3-V4
#########################

#Run R
R
#Load library
library(dada2); packageVersion("dada2")
#Define the path for the fastq files
path <- "/DATA/mikrobiom/16S_2015_DADA2/Clipped" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names

#Inspect read quality profiles
plotQualityProfile(fnFs[1:11])
plotQualityProfile(fnRs[1:11])
#For saving
library(ggplot2)
plot.quals <- plotQualityProfile(fnFs[1:11])
ggsave("qualplotFs.pdf", plot.quals, device="pdf")

plot.quals <- plotQualityProfile(fnRs[1:11])
ggsave("qualplotRs.pdf", plot.quals, device="pdf")

# Filer and trim
#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#fnFs is forward read input,filtFs is output, maxN is maximum number of ambiguous nucleotides (dada2 cannot work on amb.n.), 
#maxEE is maximal expected errors, truncQ truncates reads at the first base pair with a quality score <= 2  - 63% chance incorect
#rm.phix spike in DNA from bacteriophage calles phi X, compress and multithread speads up analyses
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,220),
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot the error rates - sanity check (we expect lower error frequency at higher quality score)
plot.err <- plotErrors(errF, nominalQ=TRUE)
ggsave("plotErrorsFs220.pdf", plot.err, device="pdf")
plot.err <- plotErrors(errR, nominalQ=TRUE)
ggsave("plotErrorsRs220.pdf", plot.err, device="pdf")

#Dereplicate sequencing (combine all identical sequencing reads into "unique sequences" with corresponding "abundance" - number of reads of that unique sequence)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inference (core sample inference algorithm - using error model developed earlier, the algorithm calc. abundance p-valuse for each unique seq,
#low p value means that is not seq error, but is actual seq)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecting the return dada-class object
dadaFs[[1]]
dadaRs[[1]]
head(getSequences(dadaFs[[1]]))

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Chimera checking and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
#Silva database https://zenodo.org/record/3731176#.XuiFLOdS9EY
#Silva files in folder 16S_2015_DADA2

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)

#Add species
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz")

#Inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#Analysis of microbiome data with phyloseq
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

save.image("16S_2015_DADA2.Rdata")

#quit R
g()