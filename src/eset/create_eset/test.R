library(testthat)
library(Biobase)

### VIASH START
meta <- list(
  resources_dir = "src/eset/create_eset/test_data",
  executable = "target/executable/eset/create_eset/create_eset"
)

### VIASH END

output <- tempfile()

out <- processx::run(meta$executable, c(
  "--pDataFile", file.path(meta$resources_dir, "pData.tsv"),
  "--fDataFile", file.path(meta$resources_dir, "fData.tsv"),
  "--mappingDir", file.path(meta$resources_dir, "mapping_dir", "AACAAGGTAC"),
  "--mappingDir", file.path(meta$resources_dir, "mapping_dir", "ACGCCTTCGT"),
  "--poolName", "foo",
  "--output", output
))
expect_equal(out$status, 0)
expect_true(file.exists(output))
result <- readRDS(output)
stopifnot(length(sampleNames(result)) == 2)
stopifnot(all(sampleNames(result) == c("foo_AACAAGGTAC", "foo_ACGCCTTCGT")))
expected_feature_names <- c(
    "ENS0001058", "ENS0000221", "ENS0001387", "ENS0000508", "ENS0001199",
    "ENS0000477", "ENS0001457", "ENS0001040", "ENS0000114", "ENS0000821",
    "ENS0001429", "ENS0001396", "ENS0000355", "ENS0000122", "ENS0000441",
    "ENS0001223", "ENS0001431", "ENS0000042", "ENS0000443", "ENS0000389",
    "ENS0001208", "ENS0001140", "ENS0000071", "ENS0001369"
)

stopifnot(length(featureNames(result)) == 24)
stopifnot(all(featureNames(result) == expected_feature_names))
expected_expressions <- matrix(
    c(0, 0,
      0, 40,
      0, 0,
      0, 0,
      1, 2,
      0, 0,
      0, 0,
      0, 0,
      2, 2,
      0, 0,
      0, 0,
      8, 2,
      0, 0,
      1, 0,
      2, 3,
      0, 0,
      0, 0,
      0, 0,
      1, 0,
      0, 0,
      16, 13,
      0, 0,
      12, 13,
      5, 2),
    ncol = 2,
    nrow = 24,
    byrow = TRUE,
)
rownames(expected_expressions) <- expected_feature_names
colnames(expected_expressions) <- c("foo_AACAAGGTAC", "foo_ACGCCTTCGT")
stopifnot(identical(exprs(result), expected_expressions))

input_f_data <- read.table(file.path(meta$resources_dir, "fData.tsv"),
                           sep = "\t", quote = "\"", comment.char = "",
                           header = TRUE)
input_f_data <- input_f_data[input_f_data$gene_id %in% expected_feature_names, ]
row.names(input_f_data) <- input_f_data$gene_id
input_f_data[] <- lapply(input_f_data, as.character)
stopifnot(identical(input_f_data, fData(result)))

# Check results filtering of barcodes with no reads
out <- processx::run(meta$executable, c(
  "--pDataFile", file.path(meta$resources_dir, "pData.tsv"),
  "--fDataFile", file.path(meta$resources_dir, "fData.tsv"),
  "--mappingDir", file.path(meta$resources_dir, "mapping_dir", "AACAAGGTAC"),
  "--mappingDir", file.path(meta$resources_dir, "mapping_dir", "EMPTY"),
  "--poolName", "bar",
  "--output", output
))
expect_equal(out$status, 0)
expect_true(file.exists(output))
result <- readRDS(output)
stopifnot(length(sampleNames(result)) == 1)
stopifnot(all(sampleNames(result) == c("bar_AACAAGGTAC")))
expected_feature_names <- c(
    "ENS0001058", "ENS0000221", "ENS0001387", "ENS0000508", "ENS0001199",
    "ENS0000477", "ENS0001457", "ENS0001040", "ENS0000114", "ENS0000821",
    "ENS0001429", "ENS0001396", "ENS0000355", "ENS0000122", "ENS0000441",
    "ENS0001223", "ENS0001431", "ENS0000042", "ENS0000443", "ENS0000389",
    "ENS0001208", "ENS0001140", "ENS0000071", "ENS0001369"
)
stopifnot(length(featureNames(result)) == 24)
stopifnot(all(featureNames(result) == expected_feature_names))
expected_expressions <- matrix(
    c(0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      2,
      0,
      0,
      8,
      0,
      1,
      2,
      0,
      0,
      0,
      1,
      0,
      16,
      0,
      12,
      5),
    ncol = 1,
    nrow = 24,
    byrow = TRUE,
)
rownames(expected_expressions) <- expected_feature_names
colnames(expected_expressions) <- c("bar_AACAAGGTAC")
stopifnot(identical(exprs(result), expected_expressions))