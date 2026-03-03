test_that("human genes successfully converted to mouse (BRCA)", {

  brca.data = selectData("wu_et_al_BRCA")
  brca.genes = sample(rownames(brca.data),300)

  converted.genes = convert_human_gene_list(brca.genes)

  # check converted genes are in data.frame or character vector format
  expect_true(is.data.frame(converted.genes) || is.vector(converted.genes))

  # check genes are in title case format or if starting with numbers, all lowercase
  check_genes_alpha = converted.genes[grepl("^[A-Za-z]",converted.genes)]
  expect_true(all(check_genes_alpha == str_to_title(check_genes_alpha)))

  check_genes_num = converted.genes[grepl("^[0-9]",converted.genes)]
  if (length(check_genes_num) > 0){
    expect_true(all(grepl("[a-z0-9]+",check_genes_num)))
  }

})
