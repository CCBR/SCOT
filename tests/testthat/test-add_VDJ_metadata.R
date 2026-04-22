# Helper: replicate the prefix logic from add_VDJ_metadata
apply_prefix_logic <- function(vdj_file) {
  prefix <- basename(vdj_file)
  if (prefix == "outs") {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }
  return(prefix)
}

test_that("prefix is empty string when basename is 'outs'", {
  expect_equal(apply_prefix_logic("/path/to/outs"), "")
  expect_equal(apply_prefix_logic("outs"), "")
})

test_that("prefix is basename with trailing underscore for non-'outs' directories", {
  expect_equal(apply_prefix_logic("/path/to/sample1"), "sample1_")
  expect_equal(apply_prefix_logic("/path/to/TCR"), "TCR_")
  expect_equal(apply_prefix_logic("/path/to/TCR_data"), "TCR_data_")
  expect_equal(apply_prefix_logic("simple_name"), "simple_name_")
})

test_that("add_VDJ_metadata errors on NULL so", {
  expect_error(add_VDJ_metadata(NULL, "some_path"))
})

test_that("add_VDJ_metadata errors on NULL vdj_file", {
  skip_if_not_installed("djvdj")
  test_so <- selectData("wu_et_al_BRCA")
  expect_error(add_VDJ_metadata(test_so, NULL))
})

test_that("add_VDJ_metadata errors when so is not a Seurat object", {
  expect_error(add_VDJ_metadata("not_a_seurat_object", "some_path"))
})
