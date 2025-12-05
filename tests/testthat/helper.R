import_pbmc <- function() {
  nFeature_RNA <- percent.mt <- NULL
  if (grepl("pbmc3k", SeuratData::InstalledData()[1]) == FALSE) {
    SeuratData::InstallData("pbmc3k")
  }

  # prepare test dataset
  pbmc3k <- suppressWarnings(SeuratData::LoadData("pbmc3k"))
  pbmc <- Seurat::UpdateSeuratObject(pbmc3k)
  set.seed(42)
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
  )
  pbmc <- Seurat::SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
  pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30)
  pbmc <- Seurat::FindClusters(pbmc, resolution = 0.8, verbose = FALSE)

  return(pbmc)
}
