
#install.packages("fs")

# 0. Recreate NAMESPACE and man/ from roxygen tags

remove.packages("OptSLDP")
.rs.restartR()

# 2. Recreate NAMESPACE and man/ from roxygen tags
source("data-raw/generate_example_data.R")

# 3. Recreate NAMESPACE and man/ from roxygen tags
devtools::document()

# 4. Now compileAttributes will work (DESCRIPTION exists, src/ exists)
Rcpp::compileAttributes()

# 5. Document again to pick up RcppExports.R
devtools::document()


# Step 1: Detach and unload the package completely
try(detach("package:OptSLDP", unload = TRUE, force = TRUE), silent = TRUE)

# Step 2: Restart R session to release the file lock
.rs.restartR()

# 6. Install
devtools::install()

# 7. Test
devtools::test()


# 8. Check
devtools::check()

# 9. Build vignettes
#options(pkgdown.internet = FALSE)

# Build everything except home, then build home separately
pkgdown::build_reference()
pkgdown::build_articles()
pkgdown::build_news()

# Build home with network disabled at the curl level
httr2_mock <- function(...) stop("no network", call. = FALSE)
pkgdown::build_home()

#unloadNamespace("OptSLDP")

#pkgdown::clean_site(force = TRUE)
pkgdown::build_site()


# 10. Build package

devtools::build()
