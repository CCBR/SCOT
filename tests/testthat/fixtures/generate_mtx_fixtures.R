# Script to generate simulated MTX fixture files for testing convert_mtx_to_h5
# Run this script from the fixtures directory to create test data
setwd(
  "/Users/bianjh/Documents/Projects/Development/SCOT/tests/testthat/fixtures"
)

library(Matrix)

# Set seed for reproducibility
set.seed(42)

# Create simulated count matrix (sparse matrix)
n_features <- 100
n_cells <- 50

# Generate sparse counts matrix
counts <- Matrix::Matrix(
  rpois(n_features * n_cells, lambda = 2),
  nrow = n_features,
  ncol = n_cells,
  sparse = TRUE
)

# Create feature names (genes)
feature_names <- data.frame(
  gene_id = paste0("ENSG_", seq_len(n_features)),
  gene_name = paste0("GENE_", seq_len(n_features))
)

# Create barcode names (cells)
barcodes <- paste0("CELL_", seq_len(n_cells))

# Set rownames and colnames
rownames(counts) <- feature_names$gene_id
colnames(counts) <- barcodes

# Write MTX file (Matrix Market format)
mtx_file <- "test_matrix.mtx.gz"
Matrix::writeMM(counts, file = sub(".gz$", "", mtx_file))

# Write barcodes file
barcodes_df <- data.frame(barcodes = barcodes)
write.table(
  barcodes_df,
  file = "test_barcodes.tsv.gz",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)

# Write features file
write.table(
  feature_names,
  file = "test_features.tsv.gz",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)
