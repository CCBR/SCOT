import_pbmc <- function() {
  if (grepl("pbmc3k", SeuratData::InstalledData()[1]) == FALSE) {
    SeuratData::InstallData("pbmc3k")
  }

  # prepare test dataset
  pbmc3k <- suppressWarnings(SeuratData::LoadData("pbmc3k"))
  pbmc <- UpdateSeuratObject(pbmc3k)
  set.seed(42)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
  )
  pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:30)
  pbmc <- FindClusters(pbmc, resolution = 0.8, verbose = FALSE)

  return(pbmc)
}
