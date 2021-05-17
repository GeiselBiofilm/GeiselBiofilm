
#install.packages("devtools")
library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") 
# change the ref argument to get other versions
#Log into server, and copy this file and RDP file into the folder with all of the samples you want to look at 
#Open R by typing R into the terminal 
#install DADA2 & DECIPHER if they are not already installed 
#Quit R
#type: nohup Rscript DADA2_DECIPHER_RDP.R &
#This will run the script - it will take a while; I recommend testing on a subset of files before doing all of them

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.11")

library("dada2"); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

#Set the filepath to your working directory
path <- getwd()
list.files(path)

#Sort forward file names
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fq.gz", full.names = TRUE))
head(fnFs)
length(fnFs)

#Sort forward file names
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fq.gz", full.names = TRUE))
head(fnRs)
length(fnRs)

#This will return a single sample name for each pair of Forward and Reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`,1)
head(sample.names)
length(sample.names)

#Visualize quality profiles of forward & reverse  reads
pdf(file = "QualF.pdf",
    width = 10,
    height = 10)
plotQualityProfile(fnFs[1:10])
dev.off()

pdf(file = "QualR.pdf",
    width = 10,
    height = 10)
plotQualityProfile(fnRs[1:10])
dev.off()

# Place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

head(filtFs)
head(filtRs)

#trim Forward reads at the 220th bp and rev reads at the 160th
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,220),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#This takes a little <10min each
#Option to increase nbases parameter
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#Check learnErrors function by plotting
#Both of these produce a warning message: Transformation introduced infinite values in continuous y-axis
pdf(file = "errF.pdf",
    width = 10,
    height = 10)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(file = "errR.pdf",
    width = 10,
    height = 10)
plotErrors(errR, nominalQ=TRUE)
dev.off()

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

#Inspecting returned data-class object
dadaFs[[1]]
dadaRs[[2]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#Construct sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
write.csv(seqtab, "seqtab.csv")

# Inspect distribution of sequence lengths
Tseq <- table(nchar(getSequences(seqtab)))
Tseq
write.csv(Tseq, "Tseqtab.csv")

#Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track.csv")

##DECIPHER ##RDP
## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## loading reference taxonomy object
load("RDP_v16-mod_March2018.RData")

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#DECIPHER, silva training set 
## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## downloading DECIPHER-formatted SILVA v138 reference
#download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DECIPHER_silva_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "DECIPHER_silva_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "DECIPHER_silva_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)


#DADA2, RDP training set 
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_16.fa.gz", multithread=TRUE)

#optional step to add species assignment 
taxa <- addSpecies(taxa, "rdp_species_assignment_16.fa.gz")

#inspect taxa results
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DADA2_RDP_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "DADA2_RDP_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "DADA2_RDP_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#DADA2, silva training set 
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

#optional step to add species assignment 
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz")

#inspect taxa results
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "DADA2_silva_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "DADA2_silva_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "DADA2_silva_ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

print("COMPLETE")
