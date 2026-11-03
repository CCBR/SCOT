SeuratData::InstallData("ifnb")

ifnb <- SeuratData::LoadData("ifnb")

ifnb_sub <- subset(
  ifnb,
  cells = colnames(ifnb)[1:200],
  features = Seurat::VariableFeatures(ifnb)[1:500]
)

saveRDS(ifnb_sub, file = test_path("fixtures", "ifnb_sub.rds"))
