
# TO DO: 
# change quality plots back to 10 samples 
# add variables to main.sh for any desired settings in pipeline.R (will need to add argument handling to this script)
# section sdript more for readability 
# collapse silva and rdp class into function?
# add quality & exploratory plots from phyloseq
# make sure ordinate plotting completes with full sample set 
# update annotation 
# add parameters for trunc length in filterandtrim

## downloading DECIPHER-formatted SILVA v138 reference
#download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

# clarifications:
# why is used trimminhg performed again


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: 
# 

# Inputs: 

# Outputs:

# Usage: 


# # Author: Owen Wilkins (Data Analytics Core), adapoted from Courtney Price & Rebecca Vals 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load required packages 
library(dada2)
library(DECIPHER)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prepare environment variables for analysis 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set directory variables 
in_dir <- getwd()
out_dir <- "./dada2/"

# read in text file containing sample names & metadata 
samples <- read.table("samples.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)

# sort samples based on name 
sample.names <- sort(samples$fastq)[10:11]

# set variables for relative path to fastq files 
fnFs <- sort(paste0("trim/", sample.names, ".R1.trim.fastq.gz"))
fnRs <- sort(paste0("trim/", sample.names, ".R2.trim.fastq.gz"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run dada2 pipeline 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# quality plots for FWD reads 
pdf(file = paste0(out_dir, "QualF.pdf"), width = 10, height = 10)
plotQualityProfile(fnFs)
dev.off()

# and for REV reads 
pdf(file = paste0(out_dir, "QualR.pdf"), width = 10, height = 10)
plotQualityProfile(fnRs)
dev.off()

# set variables for filtered fastq file names that will be created below  
filtFs <- paste0("dada2/", sample.names, "_F_filt.fastq.gz")
filtRs <- paste0("dada2/", sample.names, "_R_filt.fastq.gz")
#names(filtFs) <- sample.names
#names(filtRs) <- sample.names

# trim Forward reads at the 220th bp and rev reads at the 160th
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                      truncLen=c(220,220),
                      maxN=0, maxEE=c(2,5), 
                      truncQ=2, 
                      rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 

# learn error rates for each sample fq 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Check learnErrors output for forward reads 
pdf(file = paste0(out_dir, "errF.pdf"), width = 10, height = 10)
plotErrors(errF, nominalQ=TRUE)
dev.off()
#### Potential warning: 
# Both of these produce a warning message: Transformation introduced infinite values in continuous y-axis

# Check learnErrors output for reverse reads 
pdf(file = paste0(out_dir, "errR.pdf"), width = 10, height = 10)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# sample inference 
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table & write to file 
seqtab <- makeSequenceTable(mergers)
write.csv(seqtab, paste0(out_dir, "seqtab.csv"))

# Inspect distribution of sequence lengths & write to file 
Tseq <- table(nchar(getSequences(seqtab)))
write.csv(Tseq, paste0(out_dir, "Tseqtab.csv"))

# Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# count reads at each stage of pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# write to file 
write.csv(track, paste0(out_dir, "track.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign taxa labels 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# DECIPHER & RDP
# DECIPHER & SILVA
# DADA2 & RDP
# DADA2 & SILVA







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign taxa labels w/ DECIPHER & RDP training set
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## loading reference taxonomy object
load(paste0("references/", "RDP_v16-mod_March2018.RData"))

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
write(asv_fasta, paste0(out_dir, "ASVs.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(out_dir, "ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

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

write.table(asv_tax, paste0(out_dir.fa, "ASVs_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign taxa labels w/ DECIPHER & SILVA training set
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## loading reference taxonomy object
load(paste0("references/", "SILVA_SSU_r138_2019.RData"))

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
write(asv_fasta, paste0(out_dir, "DECIPHER_silva_ASVs.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(out_dir, "DECIPHER_silva_ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

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

write.table(asv_tax, paste0(out_dir, "DECIPHER_silva_ASVs_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)














#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign taxa labels w/ DADA2 & RDP training set 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "references/rdp_train_set_16.fa.gz", multithread=TRUE)

#optional step to add species assignment 
taxa <- addSpecies(taxa, "references/rdp_species_assignment_16.fa.gz")

#inspect taxa results
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

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
write(asv_fasta, paste0(out_dir, "DADA2_RDP_ASVs.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(out_dir, "DADA2_RDP_ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

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

write.table(asv_tax, paste0(out_dir, "DADA2_RDP_ASVs_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assign taxa labels w/ DADA2 & SILVA training set
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#DADA2, silva training set 
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "references/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

#optional step to add species assignment 
taxa <- addSpecies(taxa, "references/silva_species_assignment_v138.fa.gz")

#inspect taxa results
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

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
write(asv_fasta, paste0(out_dir, "DADA2_silva_ASVs.fa"))

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(out_dir, "DADA2_silva_ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

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

write.table(asv_tax, paste0(out_dir, "DADA2_silva_ASVs_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)










#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basic visualization w/ phyloseq package 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(Subject=samples.out)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

pdf(file = paste0(out_dir, "richness.pdf"), width = 15, height = 10)
plot_richness(ps, x="Subject", measures=c("Shannon", "Simpson"), color=1)
dev.off()

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

pdf(file = paste0(out_dir, "ordination.pdf"), width = 10, height = 10)
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
dev.off()


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

pdf(file = paste0(out_dir, "abundance_top20.pdf"), width = 15, height = 10)
plot_bar(ps.top20, x = "Subject", fill="Family")
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run checks 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print successful if all files present
# print fail if not 

print("COMPLETE")





