#install.packages("whis")
#library(testthat, warn.conflicts = FALSE)
#install.packages("processx")
#install.packages("whisker")
library(whisker)
cat(">> Running component create_report for test case \n")

out <- processx::run("./create_report", c(
  "--esetDir", "test_data",
  "--outputDir", "output"
))

expect_equal(out$status, 0)
expect_true(file.exists("output/report.Rmd"))
expect_true(file.exists("output/report.html"))
expect_true(file.exists("output/graphs"))


cat(">>  Test succesfull \n")


