library(whisker)
library(logger)

log_info("Setting temporary directory to: {meta$temp_dir}")
Sys.setenv(TMP = meta$temp_dir)
temp_folder <- tempdir(check = TRUE)
log_info("Created temporary directory {temp_folder}")

# Input esets
eset_dir <- paste(par$eset_dir, collapse = ";")
log_info("Input eset directory: {eset_dir}")


## Report: preliminary QC report
temp_rmarkdown_file <- tempfile(
  pattern = gsub(".html", "", basename(par$output_report)),
  tmpdir = temp_folder,
  fileext = ".Rmd"
)
log_info("Created temporary Rmarkdown file: {temp_rmarkdown_file}")

overall_data <- list(
  template = temp_rmarkdown_file,
  esetDir = eset_dir,
  outputDir = temp_folder
)

template <- file.path(meta$resources_dir, "template.Rmd")
log_info(paste0(
  "Rendering markdown {template} to HTML ",
  "{par$output_report}"
))
rmarkdown::render(
  file.path(meta$resources_dir, "template.Rmd"),
  output_file = basename(par$output_report),
  output_dir = dirname(par$output_report),
  params = list(
    esetDir = eset_dir,
    outputDir = temp_folder
  )
)

log_info("Done")
