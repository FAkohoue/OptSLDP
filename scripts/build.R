# Build helper script for OptSLDP

#remove.packages("fs")
#install.packages("fs")

remove.packages("OptSLDP")

# install.packages(
#   c("rlang", "fs", "cli", "glue", "lifecycle",
#     "roxygen2", "devtools", "pkgdown")
# )


# 1. Regenerate example data (fixes the 'example_sldp' vs 'example_optsldp' mismatch)
source("data-raw/generate_example_data.R")

# 2. Regenerate docs (picks up new @importFrom utils head, globalVariables)
devtools::document()

# 3. Install, test and check
devtools::install()
devtools::test()
devtools::check()


unloadNamespace("OptSLDP")

pkgdown::clean_site(force = TRUE)
pkgdown::build_site()


devtools::build()


curl::nslookup("cloud.r-project.org")

download.file("https://cloud.r-project.org", tempfile(), quiet = FALSE)

# Run from your package root in R
file.create("docs/.nojekyll")

list.files("docs", all.files = TRUE)

file.exists(".github/pkg.lock")

# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("VariantAnnotation")
