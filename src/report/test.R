library(whisker)
library(testthat)
library(R.utils)

cat(">> Creating temporary directory \n")
Sys.setenv(TMP = meta$temp_dir)
temp_folder <- tempdir(check = TRUE)

cat(">> Running component create_report for test case \n")

input_dir <- file.path(meta$resources_dir, "test_data")
stopifnot(file.exists(input_dir))


out <- processx::run(meta$executable, c(
  "--eset", file.path(meta$resources_dir, "test_data", "eset.sample_one.rds"),
  "--eset", file.path(meta$resources_dir, "test_data", "eset.sample_two.rds"),
  "--output_report", "report.html"
))

expect_equal(out$status, 0)
expect_true(file.exists("report.html"))

cat(">>  Test succesful \n")

cat(">> Running component create_report with symbolic links \n")

link_sample_1 <- file.path(temp_folder, "eset.sample_one.rds")
link_sample_2 <- file.path(temp_folder, "eset.sample_two.rds")
link_params <- file.path(temp_folder, "params.yaml")
createLink(link = link_sample_1,
           target = file.path(meta$resources_dir, "test_data", "eset.sample_one.rds"))
createLink(link = link_sample_2,
           target = file.path(meta$resources_dir, "test_data", "eset.sample_two.rds"))
createLink(link = link_params,
           target = file.path(meta$resources_dir, "test_data", "params.yaml"))

out <- processx::run(meta$executable, c(
  "--eset", link_sample_1,
  "--eset", link_sample_2,
  "--preproc_params", link_params,
  "--output_report", "report2.html"
))

expect_true(file.exists("report2.html"))