# filtering-18S
# Install if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("decontam")

library(phyloseq)
library(decontam)

setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/METABARCODING/SECUENCIACION_3/analyses2025/FILES_FEB2025/18S_feb2025")
list.files()

# Read zOTU table (rows = zOTUs, columns = samples/controls)
zOTU_table <- read.table("zOTU_table_18S_nochim.txt", 
                         header = TRUE, 
                         sep = "\t", 
                         row.names = 1, 
                         check.names = FALSE)
# Check dimensions and structure
dim(zOTU_table)

# Read metadata (rows = samples, columns = metadata fields)
metadata <- read.table("metadata_libraries_18S.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dim(metadata)
head(metadata)

# Verify all column names in zOTU_table match SampleID in metadata
all(colnames(zOTU_table) %in% metadata$SampleID)   # should return TRUE
all(metadata$SampleID %in% colnames(zOTU_table))   # should be TRUE

# Set rownames explicitly (in case not done yet)
rownames(metadata) <- metadata$SampleID

# Reorder metadata rows
metadata <- metadata[colnames(zOTU_table), ]
all(colnames(zOTU_table) == rownames(metadata))
# This should return TRUE now

##### REMOVE COUNTS LOWER THAN 4 READS----

# Set all read counts <4 to zero
zOTU_table[zOTU_table < 4] <- 0

# Verify that no counts below 4 remain (should return TRUE)
all(zOTU_table[zOTU_table > 0] >= 4)

#OPTIONAL: Remove zOTUs with zero counts across all samples
zOTU_table_h3 <- zOTU_table[rowSums(zOTU_table) > 0, ]

# Verify dimensions again
dim(zOTU_table_h3)

##### PHYLOSEQ----

# Create phyloseq object
physeq <- phyloseq(otu_table(zOTU_table_h3, taxa_are_rows = TRUE),
                   sample_data(metadata))

# Inspect object
physeq
identical(colnames(otu_table(physeq)), rownames(sample_data(physeq)))# should return TRUE

################################
### SUBSET DATA BY LIBRARIES----
### SELECT YOUR LIBRARY HERE ###
library_ID <- "LIB5"  # Change explicitly: "LIB2", "LIB3", "LIB4", or "LIB5"

### SUBSET DATA BY LIBRARY ----
physeq_lib <- subset_samples(physeq, LIB == library_ID)
physeq_lib

### DEFINE SAMPLE TYPES ----
sample_data(physeq_lib)$is.neg <- sample_data(physeq_lib)$Type == "negative_control"
sample_data(physeq_lib)$is.TSC <- sample_data(physeq_lib)$Type == "TSC"
sample_data(physeq_lib)$is.sample <- sample_data(physeq_lib)$Type == "feces"

#### NC FILTERING WITH PREVALENCE METHOD (DECONTAM)----

# Prevalence in negative controls explicitly
ps_neg <- prune_samples(sample_data(physeq_lib)$is.neg == TRUE, physeq_lib)
ps_samples <- prune_samples(sample_data(physeq_lib)$is.sample == TRUE, physeq_lib)
ps_TSC <- prune_samples(sample_data(physeq_lib)$is.TSC, physeq_lib)

prevalence_NC <- apply(otu_table(ps_neg), 1, function(x) sum(x > 0)/nsamples(ps_neg))

# Set explicit threshold (≥50% recommended)
nc_threshold <- 0.5
contaminants_NC <- names(prevalence_NC[prevalence_NC >= nc_threshold])
length(contaminants_NC)
contaminants_NC

# Calculate the MAXIMUM reads per contaminant zOTU in Negative Controls (only contaminants):
max_NC_reads <- apply(otu_table(ps_neg)[contaminants_NC, ], 1, max)

# Extract OTU tables from samples as matrix:
otu_samples <- as(otu_table(ps_samples), "matrix")

# Explicitly subtract NC maximum counts (ONLY for identified contaminants):
for (zOTU in contaminants_NC) {
  otu_samples[zOTU, ] <- otu_table(ps_samples)[zOTU, ] - max_NC_reads[zOTU]
  # Replace negative values explicitly with zero
  otu_samples[zOTU, ][otu_samples[zOTU, ] < 0] <- 0
}
# Verify explicitly no negative values remain:
all(otu_samples >= 0) # Should return TRUE

# Combine with TSCs
otu_TSC <- as(otu_table(ps_TSC), "matrix")
# Combine the corrected sample abundances with original TSC abundances explicitly
otu_samples_TSC <- cbind(otu_samples, otu_TSC)
dim(otu_samples_TSC)

# Replace original OTU table in phyloseq with the corrected abundances in samples and keeping TSCs
physeq_lib_corrected <- phyloseq(otu_table(otu_samples_TSC, taxa_are_rows = TRUE),
                                 sample_data(physeq_lib)[colnames(otu_samples_TSC), ])

# OPTIONAL remove taxa with zero total abundance after correction:
physeq_lib_corrected <- prune_taxa(taxa_sums(physeq_lib_corrected) > 0, physeq_lib_corrected)

# Verify explicitly final corrected phyloseq object
physeq_lib_corrected

########################
### TSC FILTERING BASED ON TSC RELATIVE ABUNDANCE
if (nsamples(ps_TSC) > 0) {
  
  # Extract OTU matrix from corrected phyloseq
  otu_mat <- as(otu_table(physeq_lib_corrected), "matrix")
  
  # Define TSC and sample columns
  tsc_samples <- sample_names(ps_TSC)
  real_samples <- sample_names(ps_samples)
  
  # Initialize a matrix to store filtered reads
  filtered_otu_mat <- otu_mat
  
  for (zOTU in rownames(otu_mat)) {
    
    total_reads_zOTU <- sum(otu_mat[zOTU, ])
    
    if (total_reads_zOTU > 0 && any(otu_mat[zOTU, tsc_samples] > 0)) {
      
      # Compute max rel. abundance in TSCs
      rel_abund_TSC <- otu_mat[zOTU, tsc_samples] / total_reads_zOTU
      max_rel_abund_TSC <- max(rel_abund_TSC)
      
      # Apply filter to real samples
      rel_abund_samples <- otu_mat[zOTU, real_samples] / total_reads_zOTU
      to_zero <- which(rel_abund_samples <= max_rel_abund_TSC)
      
      # Set to zero if below or equal to TSC max rel. abundance
      if (length(to_zero) > 0) {
        filtered_otu_mat[zOTU, real_samples[to_zero]] <- 0
      }
    }
  }
  
  # Create final filtered phyloseq object
  physeq_lib_final <- phyloseq(
    otu_table(filtered_otu_mat, taxa_are_rows = TRUE),
    sample_data(physeq_lib)[colnames(filtered_otu_mat), ]
  )
  
  # Remove zero taxa again (optional)
  physeq_lib_final <- prune_taxa(taxa_sums(physeq_lib_final) > 0, physeq_lib_final)
  
} else {
  # If no TSCs exist, skip this step
  message("No TSC samples available in library: ", library_ID)
  physeq_lib_final <- physeq_lib_corrected
}

# REMOVE TSCs to keep only real biological samples (feces)
physeq_lib_final <- subset_samples(physeq_lib_corrected, is.sample == TRUE)

# OPTIONAL: remove OTUs with zero total abundance again after TSC removal
physeq_lib_final <- prune_taxa(taxa_sums(physeq_lib_final) > 0, physeq_lib_final)

# Inspect final object
physeq_lib_final
  
#CHANGE TO LIBRARY NUMBER
physeq_lib5_final <-physeq_lib_final 
physeq_lib5_final

# REPEAT WITH OTHER LIBRARIES; LINE 62: SUBSET LIBRARIES

#### MERGE FILTERED LIBRARIES ----

physeq_all_filtered <- merge_phyloseq(physeq_lib1_final, physeq_lib2_final, physeq_lib3_final, physeq_lib4_final, physeq_lib5_final)
physeq_all_filtered

#### REMOVE SINGLETONS -----

# Count in how many samples each zOTU appears (>0 reads)
zOTU_sample_counts <- apply(otu_table(physeq_all_filtered), 1, function(x) sum(x > 0))

# Identify singletons (zOTUs appearing in only one sample)
singleton_zOTUs <- names(zOTU_sample_counts[zOTU_sample_counts == 1])

# Number of singletons identified explicitly
cat("Number of singleton zOTUs:", length(singleton_zOTUs), "\n")

# Remove singleton zOTUs explicitly
physeq_all_filtered_nosingletons <- prune_taxa(zOTU_sample_counts > 1, physeq_all_filtered)
physeq_all_filtered_nosingletons

#### REMOVE SAMPLES WITH NO READS ----

# Check samples with zero total reads explicitly
zero_read_samples <- sample_sums(physeq_all_filtered_nosingletons) == 0

# Number of samples explicitly identified with zero reads
num_zero_read_samples <- sum(zero_read_samples)
cat("Number of samples with zero total reads:", num_zero_read_samples, "\n")

# Names explicitly of samples with zero reads
names(zero_read_samples)[zero_read_samples]

# Remove
physeq_final <- prune_samples(sample_sums(physeq_all_filtered_nosingletons) > 0, physeq_all_filtered_nosingletons)
physeq_final

#### CONVERT TO TABLE.txt ----

physeq_final_table <- as.data.frame(otu_table(physeq_final))
dim(physeq_final_table)
write.table(physeq_final_table, "18SzOTU_tablefilteredR_final_11072025.txt", 
            sep="\t", col.names=NA, quote=FALSE)
