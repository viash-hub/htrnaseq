library(whisker)
library(logger)

log_info("Setting temporary directory to: {meta$temp_dir}")
Sys.setenv(TMP = meta$temp_dir)
temp_folder <- tempdir(check = TRUE)
log_info("Created temporary directory {temp_folder}")

template <- file.path(meta$resources_dir, "template.Rmd")

esets_normalized <- lapply(par$eset, function(eset_path) {
  return(file.path(normalizePath(dirname(eset_path)), basename(eset_path)))
})

log_info(paste0(
  "Rendering markdown {template} to HTML ",
  "{par$output_report} with esets {paste(esets_normalized, collapse = ', ')}"
))

rmarkdown::render(
  normalizePath(template),
  output_file = basename(par$output_report),
  output_dir = dirname(par$output_report),
  runtime = "static",
  intermediates_dir = par$report_dir,
  clean = TRUE,
  params = list(
    esets = esets_normalized,
    outputDir = par$report_dir
  )
)

log_info("Done")
