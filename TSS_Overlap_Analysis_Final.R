# TSS Overlap Analysis - Final
# =============================================================================
# Author: Vamsee Pillalamarri
# Email: vpillal1@jhmi.edu
# =============================================================================

# Load Libraries ---------------------------------------------------------------
# Experimental question:
# Does a gene's expression level correlate with ATAC-Seq peak height 
# at that gene's promoter (specifically the TSS)?
library(GenomicRanges)
# library(RMySQL)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Read in peak data / Create GRanges objects -----------------------------------
setwd("/users/vamsee/Desktop/Hopkins/Rotations/McCallionLab_Rotation/ATAC_plus_RNA_Seq_Analysis/individual_bams_and_peaks/")

# Read ATAC-Seq peaks for Forebrain (FB) and Midbrain (MB)
# Read in macs2 peak calls on Joint FB and MB
# Read in the 'blacklisted' regions (curated by Sarah)
FB_Joint_peaks_all <- read.table("BAMPE_Joint_FB_peaks.noheader.xls", header=T)
MB_Joint_peaks_all <- read.table("BAMPE_Joint_MB_peaks.noheader.xls", header=T)
FB_Joint_peaks_blacklisted <- read.table("Joint_FB_peaks_blacklisted.bed")
MB_Joint_peaks_blacklisted <- read.table("Joint_MB_peaks_blacklisted.bed")

# Adjust column names
colnames(FB_Joint_peaks_all) <- c("Chr","Start","End", "Length" ,"Abs_summit", 
                                  "Pileup", "X.log10.pvalue.", "fold_enrichment",
                                  "X.log10.qvalue.","name")
colnames(MB_Joint_peaks_all) <- c("Chr","Start","End", "Length" ,"Abs_summit",
                                  "Pileup", "X.log10.pvalue.", "fold_enrichment",
                                  "X.log10.qvalue.","name")
colnames(FB_Joint_peaks_blacklisted) <- c("Chr", "Start","End")
colnames(MB_Joint_peaks_blacklisted) <- c("Chr", "Start","End")

# Create GRanges Objects
FB_Joint_peaks_all <- makeGRangesFromDataFrame(FB_Joint_peaks_all, keep.extra.columns = T)
MB_Joint_peaks_all <- makeGRangesFromDataFrame(MB_Joint_peaks_all, keep.extra.columns = T)
FB_Joint_peaks_blacklisted <- makeGRangesFromDataFrame(FB_Joint_peaks_blacklisted)
MB_Joint_peaks_blacklisted <- makeGRangesFromDataFrame(MB_Joint_peaks_blacklisted)

# Adding metadata cols?
FB_Joint_peaks_blacklisted <- 
  FB_Joint_peaks_all[FB_Joint_peaks_all %in% FB_Joint_peaks_blacklisted]
MB_Joint_peaks_blacklisted <- 
  MB_Joint_peaks_all[MB_Joint_peaks_all %in% MB_Joint_peaks_blacklisted]


# Read in UCSC mm9 Ensembl Data ----------------------------------------------
setwd("/users/vamsee/Desktop/Hopkins/Rotations/McCallionLab_Rotation/ATAC_plus_RNA_Seq_Analysis/")

# TODO Add info about version # of ENSembl datatabe
# TODO Add info about the different data tables below
ENSGene_mm9_ucsc <- read.table("ENSGene_mm9_UCSC.txt", header=FALSE, sep="\t")
colnames(ENSGene_mm9_ucsc) <- c("bin","name","chrom","strand","txStart","txEnd",
                                "cdsStart","cdsEnd","exonCount","exonStarts",
                                "exonEnds","score","name2","cdsStartStat",
                                "cdsEndStat","exonFrames")

ENSGene_symbol_mm9_ucsc <- read.table("ENSGene_KnownGene_mm9_UCSC.txt", 
                                      header=FALSE, sep="\t")
colnames(ENSGene_symbol_mm9_ucsc) <- c("name","value")

ENS_Gtp_mm9_ucsc <- read.table("ENS_Gtp_mm9_UCSC.txt", header=FALSE, sep="\t")
colnames(ENS_Gtp_mm9_ucsc) <- c("gene","transcript","protein")



# Read in 'all.rpkm' RPKM RNA-Seq data -----------------------------------------
# File used = RNA-seq expression data (in RPKM) calculated for Ensembl v64 IDs
# Note: Paul made multiple "all.rpkm" files with different ID systems (ENSEMBL, REFSEQ, etc.)
setwd("/Users/vamsee/Desktop/Hopkins/Rotations/McCallionLab_Rotation/ATAC_plus_RNA_Seq_Analysis/")
all.rpkm <- read.table("All.RPKM-ensembl.v64.txt", header=TRUE)
len <- nrow(all.rpkm)
colnames(all.rpkm) <- c("ensembl_gene_id","FB.rpkm.mean","MB.rpkm.mean")

# Add gene symbol to all.rpkm
all.rpkm$mgi_symbol <- rep(".", len)
for (i in 1:nrow(all.rpkm)) {
  cat(i,'\n')
  ENST_IDs <- as.character(ENSGene_mm9_ucsc$name[ENSGene_mm9_ucsc$name2==as.character(all.rpkm$ensembl_gene_id[i])])
  all.rpkm$mgi_symbol[i] <-
    unique(as.character(ENSGene_symbol_mm9_ucsc$value[ENSGene_symbol_mm9_ucsc$name %in% ENST_IDs]))
}
# # Write updated all.rpkm to file, in the future, read it in
# # write.table(all.rpkm,'All.RPKM-ensembl.v64.withMGISymbol.txt', sep="\t", quote = F)

# Add gene biotype to all.rpkm
# NOTE: Problem exists in assigning Biotypes for genes with multiple Tx
# ENST_biotype_mm9_ucsc <- read.table("ENST_biotype_mm9_UCSC.txt", 
#                                    header=FALSE, sep="\t")
# colnames(ENST_biotype_mm9_ucsc) <- c("name","biotype")
# Do a for loop, assign biotype this way:
# ENST_biotype_mm9_ucsc$biotype[ENST_biotype_mm9_ucsc$name %in% ENST_IDs]
# NOTE: this will produce multiple biotypes per gene
# figure out a way to assign a list of biotypes (i.e. for each Tx)
# or assign an individual biotype




# Iterate through all.rpkm ----------------------------------------------------
# TODO: Add extra detail about code below
len <- nrow(all.rpkm)

# TODO: Add more info about each of the variables below
chr <- rep(".", len)
strand <- rep(".", len)
txStart <- numeric(len)
txEnd <- numeric(len)

FB_peak_hit_count <- numeric(len)
FB_peak_overlaps <- numeric(len)
FB_peak_hits <- list()
FB_peak_hit_pileup <- numeric(len)

MB_peak_hit_count <- numeric(len)
MB_peak_overlaps <- numeric(len)
MB_peak_hits <- list()
MB_peak_hit_pileup <- numeric(len)

m <- numeric(len) # number of matching records within the ENSEMBL database for this gene
ENSGene_mm9_ucsc_geneIDs <- as.character(ENSGene_mm9_ucsc$name2)
ENSGene_mm9_ucsc_chrs <- as.character(ENSGene_mm9_ucsc$chrom)
ENSGene_mm9_ucsc_strands <- as.character(ENSGene_mm9_ucsc$strand)
ENSGene_mm9_ucsc_txStarts <- as.character(ENSGene_mm9_ucsc$txStart)
ENSGene_mm9_ucsc_txEnds <- as.character(ENSGene_mm9_ucsc$txEnd)

# Iterate through all gene ids
for (i in 1:len) {
  index <- grep(all.rpkm$ensembl_gene_id[i],ENSGene_mm9_ucsc_geneIDs, 
                fixed = TRUE, value = FALSE)
  m[i] <- length(index)
  cat(i,": m length: ",m[i]," \n")
  
  if(m[i] != 0) {
    chr[i] <- unique(ENSGene_mm9_ucsc_chrs[index]) # Used to screen chrom later
    strand[i] <- unique(ENSGene_mm9_ucsc_strands[index])

    if (m[i] > 1) {
      # > 1 listed Ensembl transcripts in UCSC mm9 data
      FB_peak_hit_count_mult <- numeric(m[i])
      FB_peak_hits_mult <- numeric(m[i])
      FB_peak_hit_pileup_mult <- numeric(m[i])
      
      MB_peak_hit_count_mult <- numeric(m[i])
      MB_peak_hits_mult <- numeric(m[i])
      MB_peak_hit_pileup_mult <- numeric(m[i])
      
      # Iterate over all TSS within this gene
      for (k in 1:m[i]) {
        # Gather TSS data, create GRanges object
        TSS_bed <- ENSGene_mm9_ucsc[index[k],c("chrom","txStart","txEnd","strand","name","name2")]
        colnames(TSS_bed) <- c("Chr","Start","End","Strand","ENS_TranscriptID","ENS_GeneID")
        TSS_bed <- makeGRangesFromDataFrame(TSS_bed, keep.extra.columns = T) # TSS object created
        
        # Account for gene strand
        if (as.logical(strand(TSS_bed)=="-")) {
          # Genes on the minus strand in UCSC need to have txEnd be assigned as the TSS
          # So, if gene is on minus strand, set the TSS "region" to be txEnd-1:txEnd
          start(TSS_bed) <- end(TSS_bed)-1
        } else {
          # Genes on the positive strand in UCSC need to have txStart assigned as the TSS
          # So, if the gene is on the positive strand, set the TSS "region" to be txStart:txStart+1
          end(TSS_bed) <- start(TSS_bed)+1
        }
        
        # Overlap FB peaks with TSS
        FB_peak_hit_count_mult[k] <- countOverlaps(TSS_bed, FB_Joint_peaks_all)
        if (FB_peak_hit_count_mult[k] == 1) {
          FB_peak_hits_mult[k] <- to(findOverlaps(TSS_bed, FB_Joint_peaks_all))
          FB_peak_hit_pileup_mult[k] <- FB_Joint_peaks_all$Pileup[FB_peak_hits_mult[k]]
        } else {
          cat("Too many peak hits within a multiple transcript Gene Model! Or NO HITS!\n")
        }
        
        # Overlap MB peaks with TSS
        MB_peak_hit_count_mult[k] <- countOverlaps(TSS_bed, MB_Joint_peaks_all)
        if (MB_peak_hit_count_mult[k] == 1) {
          MB_peak_hits_mult[k] <- to(findOverlaps(TSS_bed, MB_Joint_peaks_all))
          MB_peak_hit_pileup_mult[k] <- MB_Joint_peaks_all$Pileup[MB_peak_hits_mult[k]]
        } else {
          cat("Too many peak hits within a multiple transcript Gene Model! Or NO HITS!\n")
        }
      }
      
      # Assign unique peak overlaps only for this gene, if there are hits for any of the TSSs
      if (any(FB_peak_hit_count_mult !=0)) {
        FB_peak_hits[[i]] <- unique(FB_peak_hits_mult[FB_peak_hits_mult!=0]) # Only record unique peak overlaps amongst all TSSs
        FB_peak_hit_pileup[i] <- sum(unique(FB_peak_hit_pileup_mult[FB_peak_hits_mult!=0])) # Sum up pileup values for unique peak overlaps amongst all TSSs
        FB_peak_hit_count[i] <- length(unique(FB_peak_hits_mult[FB_peak_hits_mult!=0])) # Count up num unique peaks across all TSSs
      }
      
      if (any(MB_peak_hit_count_mult !=0)) {
        MB_peak_hits[[i]] <- unique(MB_peak_hits_mult[MB_peak_hits_mult!=0]) # Only record unique peak overlaps amongst all TSSs
        MB_peak_hit_pileup[i] <- sum(unique(MB_peak_hit_pileup_mult[MB_peak_hits_mult!=0])) # Sum up pileup values for unique peak overlaps amongst all TSSs
        MB_peak_hit_count[i] <- length(unique(MB_peak_hits_mult[MB_peak_hits_mult!=0])) # Count up num unique peaks across all TSSs
      }
      
    } else {
      # Only one isoform (thus 1 listed Ensembl transcript in UCSC mm9)
      # Gather TSS data from UCSC Ensembl data, and create a GRanges object
      # Intersect GRanges object with peak data (MB and MB separately)
      
      # Gather TSS data, create GRanges object
      TSS_bed <- ENSGene_mm9_ucsc[index,c("chrom","txStart","txEnd","strand","name","name2")]
      colnames(TSS_bed) <- c("Chr","Start","End","Strand","ENS_TranscriptID","ENS_GeneID")
      TSS_bed <- makeGRangesFromDataFrame(TSS_bed, keep.extra.columns = T)
      
      # Account for gene strand
      if (as.logical(strand(TSS_bed)=="-")) {
        # Genes on the minus strand in UCSC need to have txEnd be assigned as the TSS
        # So, if gene is on minus strand, set the TSS "region" to be (txEnd-1):txEnd
        start(TSS_bed) <- end(TSS_bed)-1
      } else {
        # Genes on the positive strand in UCSC need to have txStart assigned as the TSS
        # So, if the gene is on the positive strand, set the TSS "region" to be txStart:(txStart+1)
        end(TSS_bed) <- start(TSS_bed)+1
      }
      
      # Overlap FB peaks with TSS
      FB_peak_hit_count[i] <- countOverlaps(TSS_bed, FB_Joint_peaks_all) # Actual INTERSECTION happens here
      if (FB_peak_hit_count[i] > 0) {
        if (FB_peak_hit_count[i] == 1) {
          FB_peak_hits[[i]] <- to(findOverlaps(TSS_bed, FB_Joint_peaks_all))
          FB_peak_hit_pileup[i] <- FB_Joint_peaks_all$Pileup[FB_peak_hits[[i]]]
          FB_peak_hit_count[i] <- 1
        } else {
          cat("Too many hits for SINGLE Transcript!\n")
        }
      }
      
      # Overlap MB peaks with TSS
      MB_peak_hit_count[i] <- countOverlaps(TSS_bed, MB_Joint_peaks_all)
      if (MB_peak_hit_count[i] > 0) {
        if (MB_peak_hit_count[i] == 1) {
          MB_peak_hits[[i]] <- to(findOverlaps(TSS_bed, MB_Joint_peaks_all))
          MB_peak_hit_pileup[i] <- MB_Joint_peaks_all$Pileup[MB_peak_hits[[i]]]
          MB_peak_hit_count[i] <- 1
        } else {
          cat("Too many hits for SINGLE Transcript!\n")
        }
      }
    }
  }
}

# Create a results data frame ----
# Data frame columns: Gene (ENS..GID), MGI Symbol, Gene Biotype, MB RPKM Mean, MB Pileup@TSS, ...
#           MB RPKM Mean, MB Pileup@TSS, Tx Subtype (Single / Multiple / NA), 
#           Match (m) in UCSC Ensembl Db

# Annotate 'multiple' or 'single' transcripts for genes
tx_subtype <- character(len)
for (i in 1:len) {
  if (m[i]==0) {
    tx_subtype[i] <- "NA" # gene was not matched in UCSC Ensembl data
  } else if (m[i]==1) {
    tx_subtype[i] <- "single"
  } else {
    tx_subtype[i] <-"multiple"
  }
}
num_tx_subtypes <- c(length(which(tx_subtype=="single")),
                     length(which(tx_subtype=="multiple")),
                     length(which(tx_subtype=="NA")))

df <- data.frame(all.rpkm$ensembl_gene_id,all.rpkm$mgi_symbol, chr, strand,
                 all.rpkm$FB.rpkm.mean,FB_peak_hit_pileup,
                 all.rpkm$MB.rpkm.mean,MB_peak_hit_pileup, tx_subtype, m)
# Plot Set (FB, MB) by ggplot --------------------------------------------------
#   Plot log pileups (x) vs log mean rpkms (y) (ALL GENES)
#   Plot x vs y where 'x' is pileup values for genes (one TSS and multiple TSSs)
#   and 'y' is the pileup value for the peak that matched that TSS

df <- df[df$m>=1 & df$chr !="chrM", ] # Ignore chrM
df <- df[grep("random", df$chr, invert = TRUE),] # Ignore 'random' chr contigs
df <- df[grep("Un", df$chr, invert = TRUE), ] # Ignore 'Unplaced' chr contigs
# df <- df[!(df$all.rpkm.FB.rpkm.mean == 0 | df$FB_peak_hit_pileup == 0),] # Ignore genes with zero-RPKM OR zero-pileups

# TODO ADD Pseudocount / make it clear that "0" RPKMS / Pileups are actually 0s not other values

FB <- ggplot(df[which(df$m>=1 & df$chr!="chrM"),], aes(x=log2(df$FB_peak_hit_pileup[which(df$m>=1 & df$chr!="chrM")]), y=log2(df$all.rpkm.FB.rpkm.mean[which(df$m>=1 & df$chr!="chrM")]))) +
  geom_point(alpha=0.5,aes(colour = factor(df$tx_subtype[which(df$m>=1 & df$chr!="chrM")]))) +
  labs(x="log2 Pileup",y="log2 Mean RPKM", 
       title="Log Transformed Pileups vs Mean RPKMs at TSSs \n(n=37,681 Genes) for Forebrain") +
  scale_color_brewer(type="qual", palette="Set1", name="# Transcripts") +
  geom_vline(xintercept = 10, linetype="longdash",colour="red") +
  geom_hline(yintercept = 10, linetype="longdash",colour="red") + 
  scale_y_continuous(breaks=-10:10) + 
  #coord_cartesian(ylim = c(-10, 12), xlim=c(3,10))
  stat_smooth(method="lm", col = "yellow", fullrange = T) 

MB <- ggplot(df[which(df$m>=1 & df$chr!="chrM"),], aes(x=log2(df$MB_peak_hit_pileup[which(df$m>=1 & df$chr!="chrM")]), y=log2(df$all.rpkm.MB.rpkm.mean[which(df$m>=1 & df$chr!="chrM")]))) +
  geom_point(alpha=0.5,aes(colour = factor(df$tx_subtype[which(df$m>=1 & df$chr!="chrM")]))) +
  labs(x="log2 Pileup",y="log2 Mean RPKM", 
       title="Log Transformed Pileups vs Mean RPKMs at TSSs \n(n=37,681 Genes) for Midbrain") +
  scale_color_brewer(type="qual", palette="Set1", name="# Transcripts") +
  geom_vline(xintercept = 10, linetype="longdash",colour="red") +
  geom_hline(yintercept = 10, linetype="longdash",colour="red") + 
  scale_y_continuous(breaks=-10:10) +
  #coord_cartesian(ylim = c(-10, 12), xlim=c(3,10))
  stat_smooth(method="lm", col = "yellow", fullrange = T)

grid.arrange(MB, FB, ncol=2)
# plot_grid(MB, FB, labels = c("A", "B"), align = "h") # using cowplot


# Correlation analyses ---------------------------------------------------------
# Correlation between TSS peaks and expression levels
# Only for genes that have non-zero RPKM
df <- df[df$m>=1 & df$chr !="chrM", ] # Ignore chrM
df <- df[grep("random", df$chr, invert = TRUE),] # Ignore 'random' chr contigs
df <- df[grep("Un", df$chr, invert = TRUE), ] # Ignore 'Unplaced' chr contigs

# r^2 for FB
fb_table <- df[!(df$all.rpkm.FB.rpkm.mean == 0 | df$FB_peak_hit_pileup == 0),] # Ignore genes with zero-RPKM OR zero-pileups
pileups <- fb_table$FB_peak_hit_pileup
fpkms <- fb_table$all.rpkm.FB.rpkm.mean
z2 <- data.frame(pileups,fpkms)
r.squared <- summary(lm(log2(pileups)~log2(fpkms), data=z2))$r.squared
r.squared

# r^2 for MB
mb_table <- df[!(df$all.rpkm.MB.rpkm.mean == 0 | df$MB_peak_hit_pileup == 0),] # Ignore genes with zero-RPKM OR zero-pileups
pileups <- mb_table$MB_peak_hit_pileup
fpkms <- mb_table$all.rpkm.MB.rpkm.mean
z2 <- data.frame(pileups,fpkms)
r.squared <- summary(lm(log2(pileups)~log2(fpkms), data=z2))$r.squared
r.squared

# Linear regression of rpkms ~ pileups (OLD CODE)
#   ...with pseudo count of 1 for 0 values prior to log-transformation
# a <- df$all.rpkm.FB.rpkm.mean[which(df$m>=1)] # RPKM
# a[which(a==0)] <- 1 # This is wrong...add pseudocount to *all* values not just the ones with 0 values in them
# b <- df$FB_peak_hit_pileup[which(df$m>=1)] # Pileups
# b[which(b==0)] <- 1
# a <- log2(a)
# b <- log2(b)
# 
# z2 <- data.frame(a,b)
# 
# fit <- lm(a,b, data = z2) # Linear regression of RPKM ~ Pileups
# summary(fit2)