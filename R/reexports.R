#' walrus operator
#' @importFrom rlang :=
#' @export
rlang::`:=`

#' rlang data pronoun
#' @importFrom rlang .data
#' @export
rlang::.data


# Suppress R CMD check note 'All declared Imports should be used'.
# These packages are used by Seurat, but are only listed in Seurat's Suggests field.
# https://community.rstudio.com/t/how-should-a-meta-package-handle-this-note-all-declared-imports-should-be-used/23400/3
#' @importFrom harmony RunHarmony
#' @importFrom rliger createLiger
NULL
