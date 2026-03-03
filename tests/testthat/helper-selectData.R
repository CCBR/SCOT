selectData <- function(dataset) {

  if (dataset == "wu_et_al_BRCA"){

    print("selected wu_et_al_BRCA dataset")
    input.object <- readRDS(test_path("fixtures/",
                                      "wu_et_al_BRCA.rds"))

  }

}
