# BRCA Dataset Download Guide

This document describes how to download the BRCA single-cell RNA-seq dataset from the publication:

> Integrated single-cell and spatial transcriptomics reveals breast cancer heterogeneity
> *Nature Genetics* (2021) | [DOI: 10.1038/s41588-021-00911-1](https://www.nature.com/articles/s41588-021-00911-1)

## Data Availability

The BRCA dataset is available through multiple sources:

### 1. Broad Institute Single Cell Portal (Recommended for Interactive Exploration)

**URL:** https://singlecell.broadinstitute.org/single_cell/study/SCP1039

**Features:**
- In-browser data exploration
- Interactive visualization
- Direct download of processed data
- No registration required

**Data Type:** Processed scRNA-seq data

### 2. Gene Expression Omnibus (GEO)

**Accession:** GSE176078
**URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

**Features:**
- Curated data repository
- Multiple download formats
- Supplementary files and metadata
- NCBI account recommended (optional)

**Data Type:** Processed scRNA-seq data

### 3. European Genome-phenome Archive (EGA)

**Accession:** EGAS00001005173
**URL:** https://ega-archive.org/

**Features:**
- European Bioinformatics Institute & Centre for Genomic Regulation
- Controlled access for raw data
- May require registration and data access approval

**Data Type:** Raw scRNA-seq data

## Processing Information

The downloaded dataset was processed according to the **Seurat SCTransform and Merge workflow**:

- **Reference:** [Seurat PBMC 3K Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- **Key steps:**
  - SCTransform normalization
  - Data merge from multiple samples
  - Quality control filtering
  - Dimensionality reduction (PCA, UMAP)
  - Clustering

## Quick Start: Downloading from Broad Institute Portal

1. Navigate to https://singlecell.broadinstitute.org/single_cell/study/SCP1039
2. Browse the interactive visualizations
3. Click the **Download** button to access processed data
4. Select your desired file format (HDF5, CSV, Seurat object, etc.)
5. Save to your local machine

## Loading in R

Once downloaded, you can load the data in R:

```r
library(Seurat)

# If downloaded as a Seurat object (.rds)
seurat_obj <- readRDS("path/to/BRCA_data.rds")

# If downloaded as HDF5
library(Seurat)
seurat_obj <- Read10X_h5("path/to/BRCA_data.h5")

# If downloaded as expression matrix + metadata
expr_matrix <- read.csv("path/to/expression_matrix.csv", row.names = 1)
metadata <- read.csv("path/to/metadata.csv", row.names = 1)
seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = metadata)
```

## Data Access Considerations

- **Processed data** (Broad Portal, GEO): Generally open access, no approval needed
- **Raw data** (EGA): May require registration and data access approval (check study-specific requirements)
- **Computational requirements:** The full dataset is large; consider downsampling for initial exploration
