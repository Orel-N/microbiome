#################################
#dada2 16S Microbiome 2020
#################################

#Run R
R
#Load library
library(dada2); packageVersion("dada2")
#Define the path for the fastq files
path <- "/DATA/mikrobiom/16S_2020_DADA2/Clipped" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names

library(ggplot2)
#Inspect read quality profiles
plotQualityProfile(fnFs[1:40])
ggsave("qualplotFs_2020.pdf", last_plot())
plotQualityProfile(fnRs[1:40])
ggsave("qualplotRs_2020.pdf", last_plot())



save.image("16S_2015_DADA2.Rdata")

#quit R
q()