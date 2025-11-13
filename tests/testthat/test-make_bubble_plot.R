test_that("bubble plot is successful for (BRCA)", {
  
  brca.data = selectData("wu_et_al_BRCA")
  brca.genes = sample(rownames(brca.data),5)
  
  bubble.plot = make_bubble_plot(so = brca.data,
                                 features = brca.genes,
                                 ident = "subtype")
  
  skip_on_ci()
  expect_snapshot_file(
    .make_bubble_plot(bubble.plot),
    "brca_bubbleplot.png"
  )
  
  expected_elements = c("ggplot2::ggplot","ggplot","ggplot2::gg","S7_object","gg")
  expect_setequal(class(bubble.plot), expected_elements)
  
})