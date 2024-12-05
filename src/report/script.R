library(whisker)

outputDir <- par$outputDir
dir.create(outputDir)
esetDir <- paste(par$esetDir, collapse=";")

outputFile <- "report.html"

## Report: preliminary QC report 
reportTemplate <- "template.Rmd"
reportTemplateText <- readLines(reportTemplate)

overallData <- list(
  template = gsub(".html", ".Rmd",  outputFile),
  esetDir = esetDir,
  outputDir = outputDir)

writeLines(whisker.render(reportTemplateText, overallData), gsub(".html", ".Rmd", outputFile))

rmarkdown::render(gsub(".html", ".Rmd",  outputFile), output_file = outputFile, output_dir = outputDir)

#file.rename(gsub(".html", ".Rmd",  outputFile), gsub(".html", ".Rmd", file.path(reportDir, outputFile)))
#file.rename("graphs",  file.path(reportDir, "graphs"))
