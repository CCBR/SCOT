# Script to downsample BRCA dataset for unit testing
#
# This script reduces the BRCA dataset to a minimal size suitable for unit testing
# while preserving the structure needed for comprehensive testing.
# Target: ~1MB on disk

library(SeuratObject)
library(Seurat)

# Load the current fixture
brca <- readRDS(
  test_path(
    "fixtures",
    "BRCA_Combine_and_Renormalize_SO_downsample.rds"
  )
)

brca <- subset(
  brca,
  cells = sample(colnames(brca), 250),
  features = VariableFeatures(brca)[1:250]
)

# Remove SCT misc (vst.out parameters) ~ 711kb
# keeps the SCT assay data
if ("SCT" %in% names(brca@assays)) {
  brca@assays$SCT@misc <- list()
}

# compress
saveRDS(
  brca,
  "tests/testthat/fixtures/BRCA_Combine_and_Renormalize_SO_downsample_1mb.rds",
  compress = "xz"
)
