# Create simulated Seurat object for testing

# Create a minimal count matrix with genes and cells
set.seed(42)
n_features <- 100
n_cells <- 150

# Generate count matrix with realistic values
counts <- matrix(
  rpois(n_features * n_cells, lambda = 5),
  nrow = n_features,
  ncol = n_cells
)
rownames(counts) <- paste0("GENE_", seq_len(n_features))
colnames(counts) <- paste0("CELL_", seq_len(n_cells))

# Create Seurat object
test_seurat <- SeuratObject::CreateSeuratObject(counts = counts)

# Add metadata for batch correction (no NA/NaN values)
test_seurat@meta.data$batch <- factor(
  rep(c("batch1", "batch2"), length.out = ncol(test_seurat))
)
test_seurat@meta.data$nFeature_RNA <- colSums(counts > 0)
test_seurat@meta.data$nCount_RNA <- colSums(counts)
test_seurat@meta.data$percent_mt <- runif(ncol(test_seurat), min = 0, max = 10)

# Set identities to batch for harmony integration
SeuratObject::Idents(test_seurat) <- test_seurat@meta.data$batch

test_that("run_batch_correction returns a Seurat object", {
  corrected_result <- run_batch_correction(
    so_in = test_seurat,
    npcs = 3,
    resolution_list = c(0.5, 1.0),
    species = "hg38",
    method_in = "RPCAIntegration"
  )

  expect_s4_class(corrected_result, "Seurat")
  expect_true(ncol(corrected_result) > 0)
  expect_true(nrow(corrected_result) > 0)
})

test_that("run_batch_correction errors on invalid species", {
  expect_error(
    run_batch_correction(
      so_in = test_seurat,
      species = "invalid_species",
      method_in = "LIGER"
    ),
    regexp = "species|invalid|not supported"
  )
})

test_that("run_batch_correction errors on invalid method", {
  expect_error(
    run_batch_correction(
      so_in = test_seurat,
      species = "hg38",
      method_in = "invalid_method"
    ),
    regexp = "method|invalid|not supported"
  )
})

test_that("run_batch_correction errors on NULL Seurat object", {
  expect_error(
    run_batch_correction(
      so_in = NULL,
      species = "hg38",
      method_in = "RPCAIntegration"
    ),
    regexp = "Seurat|NULL|NULL object"
  )
})
