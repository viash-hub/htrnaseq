library(whisker)
library(testthat)
cat(">> Running component create_report for test case \n")

input_dir <- file.path(meta$resources_dir, "test_data")
stopifnot(file.exists(input_dir))

out <- processx::run(meta$executable, c(
  "--eset_dir", file.path(meta$resources_dir, "test_data"),
  "--output_report", "output.html"
))

expect_equal(out$status, 0)
expect_true(file.exists("output.html"))

cat(">>  Test succesful \n")
