load_plain_pbmc <- function() {
  if (grepl("pbmc3k", SeuratData::InstalledData()[1]) == FALSE) {
    SeuratData::InstallData("pbmc3k")
  }

  # prepare test dataset
  pbmc3k <- suppressWarnings(SeuratData::LoadData("pbmc3k"))
  pbmc <- Seurat::UpdateSeuratObject(pbmc3k)
  return(pbmc)
}

load_processed_pbmc <- function() {
  nFeature_RNA <- percent.mt <- NULL
  pbmc <- load_plain_pbmc()
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(
    pbmc,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
  )
  pbmc <- suppressWarnings(Seurat::SCTransform(
    pbmc,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  ))
  pbmc <- Seurat::RunPCA(
    pbmc,
    features = Seurat::VariableFeatures(object = pbmc)
  )
  pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30)
  pbmc <- Seurat::FindClusters(pbmc, resolution = 0.8, verbose = FALSE)

  return(pbmc)
}

load_fixture_data <- function(dataset) {
  if (dataset == "wu_et_al_BRCA") {
    return(readRDS(testthat::test_path(
      "fixtures/",
      "BRCA_Combine_and_Renormalize_SO_downsample_1mb.rds"
    )))
  }
}

#' Returns a random counts matrix.
#' Source: https://github.com/satijalab/seurat/blob/main/tests/testthat/test_dimensional_reduction.R
get_random_counts <- function() {
  # Populate a 100 by 100 matrix with random integers from 1 to 50.
  counts <- matrix(
    data = sample(c(1:50), size = 1e4, replace = TRUE),
    ncol = 100,
    nrow = 100
  )

  # Assign column and row names to the matrix to label cells and genes.
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  row.names(counts) <- paste0("gene", seq_len(nrow(counts)))

  # Convert `counts` to a `dgCMatrix`.
  counts_sparse <- Seurat::as.sparse(counts)

  return(counts_sparse)
}

#' Returns a `Seurat` instance containing the specified `assay_version` and
#' populated with `counts` which is also preprocessed (normalized + scaled).
#'
#' Source: https://github.com/satijalab/seurat/blob/main/tests/testthat/test_dimensional_reduction.R
get_test_data <- function(
  counts = get_random_counts(),
  assay_version = getOption("Seurat.object.assay.version")
) {
  # Use the `assay_version` param to choose the correct assay builder.
  create_assay <- switch(
    assay_version,
    v5 = SeuratObject::CreateAssay5Object,
    stop("`assay_version` should be 'v5'")
  )
  # And then instantiate the specified assay type.
  assay <- create_assay(counts)

  # Instantiate a `Seurat` instance using the default assay name.
  test_data <- Seurat::CreateSeuratObject(assay)

  # Normalize, and then scale the input data.
  test_data <- Seurat::NormalizeData(test_data, verbose = FALSE)
  test_data <- Seurat::ScaleData(test_data, verbose = FALSE)
  test_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    test_data,
    pattern = "^MT-"
  )

  test_data <- suppressWarnings(Seurat::SCTransform(
    test_data,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  ))
  test_data <- Seurat::RunPCA(
    test_data,
    features = Seurat::VariableFeatures(object = test_data)
  )
  test_data <- Seurat::FindNeighbors(test_data, dims = 1:30)
  test_data <- Seurat::FindClusters(
    test_data,
    resolution = 0.8,
    verbose = FALSE
  )

  return(test_data)
}

load_mouse_bm <- function() {
  set.seed(42)
  init_data <- scRNAseq::fetchDataset(
    name = "grun-bone_marrow-2016",
    version = "2023-12-14"
  )
  #Conversion from gene__chrN annotation
  gene_names <- make.names(gsub("_.*", "", rownames(init_data)), unique = TRUE)
  SummarizedExperiment::assay(init_data, "counts") <-
    as(SummarizedExperiment::assay(init_data, "counts"), "dgCMatrix")

  seurat_object <- Seurat::CreateSeuratObject(
    counts = init_data@assays@data$counts
  )

  rownames(seurat_object) <- gene_names

  return(seurat_object)
}

load_mouse_processed_bm <- function() {
  mm_bm <- load_mouse_bm()
  nFeature_RNA <- percent.mt <- NULL
  mm_bm[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    mm_bm,
    pattern = "^mt-"
  )
  mm_bm <- subset(
    mm_bm,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
  )
  mm_bm <- suppressWarnings(Seurat::SCTransform(
    mm_bm,
    verbose = FALSE
  ))
  mm_bm <- Seurat::RunPCA(
    mm_bm,
    features = Seurat::VariableFeatures(object = mm_bm)
  )
  mm_bm <- Seurat::FindNeighbors(mm_bm, dims = 1:30)
  mm_bm <- Seurat::FindClusters(mm_bm, resolution = 0.8, verbose = FALSE)

  return(mm_bm)
}
