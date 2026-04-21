selectData <- function(dataset) {
  if (dataset == "wu_et_al_BRCA") {
    return(readRDS(testthat::test_path(
      "fixtures/",
      "BRCA_Combine_and_Renormalize_SO_downsample.rds"
    )))
  }
}
