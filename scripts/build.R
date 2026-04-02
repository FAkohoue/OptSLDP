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


options(pkgdown.internet = FALSE)
pkgdown::build_site()
