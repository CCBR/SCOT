# Script to download and downsample BRCA dataset for unit testing
#
# This script downloads the BRCA single-cell RNA-seq dataset from GEO (GSE176078)
# and downsamples it to a minimal size suitable for unit testing while preserving
# the structure needed for comprehensive testing.
# Target: ~1MB on disk
#
# Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
# Publication: Wu et al. (2021), Nature Genetics

library(Seurat)
library(GEOquery)

# Create temporary directory for downloads
temp_dir <- tempdir()
geo_dir <- file.path(temp_dir, "GSE176078")
dir.create(geo_dir, showWarnings = FALSE)

cat("Downloading BRCA dataset from GEO (GSE176078)...\n")

# Download supplementary files from GEO
# This downloads the main scRNA-seq tar.gz file (~533 MB)
# Note: This may take several minutes depending on connection speed
geo_data <- getGEOSuppFiles(
  "GSE176078",
  baseDir = geo_dir,
  makeDirectory = FALSE
)

# Extract the tar.gz file
tar_file <- file.path(geo_dir, "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")

if (file.exists(tar_file)) {
  cat("Extracting compressed archive...\n")
  untar(tar_file, exdir = geo_dir)
} else {
  stop("Failed to download tar.gz file from GEO")
}

# Find the Seurat object or count matrices
# The archive typically contains:
# - A Seurat RDS object, or
# - Count matrices (barcodes.tsv, features.tsv, matrix.mtx)
rds_files <- list.files(
  geo_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(rds_files) > 0) {
  brca_full_path <- rds_files[1]
  cat("Found Seurat object:", basename(brca_full_path), "\n")
} else {
  stop(
    "Could not find Seurat object (.rds) in extracted archive.\n",
    "Files in directory:\n",
    paste(list.files(geo_dir, recursive = TRUE), collapse = "\n")
  )
}

# Load the full dataset
cat("Loading full BRCA dataset...\n")
brca <- readRDS(brca_full_path)

cat("Downsampling for testing...\n")

# Downsample for testing:
# - Subset to 250 cells (from potentially 100k+)
# - Subset to 250 variable features (from potentially 2k-3k)
brca <- subset(
  brca,
  cells = sample(colnames(brca), 250),
  features = VariableFeatures(brca)[1:250]
)

# Remove SCT misc (vst.out parameters) to reduce file size (~711kb)
# This keeps the SCT assay data intact
if ("SCT" %in% names(brca@assays)) {
  brca@assays$SCT@misc <- list()
}

# Save the downsampled dataset
# Using xz compression to minimize disk footprint
fixture_path <- test_path("fixtures", "BRCA_downsample.rds")
saveRDS(brca, fixture_path, compress = "xz")

cat("✓ BRCA test fixture created successfully\n")
cat("  Location:", fixture_path, "\n")
cat("  File size:", format(file.size(fixture_path), units = "auto"), "\n")

# Clean up temporary files
unlink(geo_dir, recursive = TRUE)
