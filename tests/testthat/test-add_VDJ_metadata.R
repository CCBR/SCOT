test_so <- load_fixture_data("wu_et_al_BRCA")

test_that("add_VDJ_metadata errors when so is not a Seurat object", {
  expect_error(
    add_VDJ_metadata("not_a_seurat_object", "some_path"),
    regexp = "does not exist"
  )
})

test_that("add_VDJ_metadata errors on NULL vdj_file", {
  expect_error(
    add_VDJ_metadata(test_so, NULL),
    regexp = "character vector argument expected|basename"
  )
})
