---
name: "R Standards"
description: "Coding conventions for R files"
applyTo: "**/*.R"
---

# R coding standards

- R scripts must include function and class docstrings via roxygen2.
- CLIs must be defined using the `argparse` package.
- CLIs must support `--help` and document required/optional arguments.
- R code should pass `lintr` and `air`.
- Tests should be written with `testthat`.
- Packages should pass `devtools::check()`.
- R code should adhere to the tidyverse style guide. https://style.tidyverse.org/
- Object names (functions, variables) should follow the snake_case style.
- Function names should start with a verb.
- Use the native pipe operator `|>`
- If an external package is required, prefer packages from the tidyverse, r-lib, or other posit-affiliated organizations.
- Prefer using tidyverse functions rather than base R functions where possible.
- Only include one return statement at the end of a function, if a return statement is used at all. Explicit returns are preferred but not required for R functions.
