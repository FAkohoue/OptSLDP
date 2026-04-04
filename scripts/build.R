# Build helper script for OptSLDP

install.packages("rsvg")
rsvg::rsvg_png(
  "man/figures/logo.svg",
  "man/figures/logo.png"
)

library(magick)

# Load your logo
img <- image_read("man/figures/logo.png")

# Create favicon sizes
sizes <- c(16, 32, 48, 64, 180, 192, 512)

dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)

for (s in sizes) {
  resized <- image_resize(img, paste0(s, "x", s))
  image_write(resized, paste0("pkgdown/assets/favicon-", s, ".png"))
}

# Create favicon.ico (multi-size)
ico <- image_resize(img, "64x64")
image_write(ico, "pkgdown/assets/favicon.ico")

writeLines(
  c(
    '<link rel="icon" type="image/svg+xml" href="logo.svg">',
    '<link rel="icon" type="image/png" sizes="32x32" href="favicon-32.png">',
    '<link rel="apple-touch-icon" href="favicon-180.png">'
  ),
  "pkgdown/assets/favicon.html"
)


# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")


#remove.packages("fs")
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
