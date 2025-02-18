library(Biobase)
library(data.table)
library(nlcv)
library(Matrix)
library(Seurat)

### VIASH START
par <- list(
  pDataFile = "src/eset/create_eset/test_data/pData.tsv",
  fDataFile = "src/eset/create_eset/test_data/fData.tsv",
  studyType = "Standard",
  mappingDir = c("src/eset/create_eset/test_data/mapping_dir/AACAAGGTAC",
                 "src/eset/create_eset/test_data/mapping_dir/ACGCCTTCGT"),
  output = "eset.rds",
  poolName = "Foo"
)
### VIASH END


Read10X <- function(data_dir = NULL, gene_column = 2, unique_features = TRUE) {
  full.data <- list()
  for (i in seq_along(along.with = data_dir)) {
    run <- data_dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "features.tsv")
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data_dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data_dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                              yes = gene.loc, no = features.loc), header = FALSE, 
                                stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene_column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene_column]))
      replacement.column <- ifelse(test = gene_column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene_column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique_features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene_column) {
        stop(paste0("gene_column was set to ", gene_column,
                    " but feature.tsv.gz (or genes.tsv) only has ",
                    fcols, " columns.", " Try setting the gene_column ",
                    "argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene_column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message(paste0("10X data contains more than one type and is ",
                       "being returned as a list containing matrices ",
                       "of each type."))
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    } else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "CsparseMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

match_features <- function(exprs_matrix, fdata) {

  identical_features <- all(rownames(exprs_matrix) == rownames(fdata))

  if (nrow(exprs_matrix) != nrow(fdata) || !identical_features) {
    message(paste0("Features in 'fData' and expression matrix differ. ",
                   "Only matching features are returned."))
  }

  features <- intersect(rownames(exprs_matrix), rownames(fdata))
  exprs_matrix <- exprs_matrix[which(rownames(exprs_matrix) %in% features), ]
  fdata <- fdata[which(rownames(fdata) %in% features), ]

  fdata[, seq_len(ncol(fdata))] <- lapply(fdata[, seq_len(ncol(fdata)), drop = FALSE], as.character)
  # order features in exprs mat according to fdata
  exprs_matrix <- exprs_matrix[match(rownames(fdata), rownames(exprs_matrix)), ]

  list(exprs_matrix = exprs_matrix, fdata = fdata)

}


create_pdata <- function(sample_file, pool_name, barcodes) {
  cols_to_remove <- c("SampleFileName", "Output", "Measure", "Strandedness")
  pData <- sample_file[, !colnames(sample_file) %in% cols_to_remove,
                       drop = FALSE]
  rownames(pData) <- lapply(sample_file$WellBC,
                            \(x) paste(pool_name, x, sep = "_"))
  # pData[, ] <- lapply(pData, as.factor)
  pData$PoolName <- pool_name
  pData <- pData[match(barcodes, pData$WellBC), ]
  return(pData)
}

check_sample_file <- function(mapping_dir, sample_file){

  message("Checking sample annotation:")

  requireNamespace("tools")
  mapping_dir <- unlist(lapply(mapping_dir, function(x) {
    if (!dir.exists(x)) {
      stop(sprintf(paste0("Could not find directory ",
                          "provided in 'mappingDir' argument (%s)."), x))
    }
    tools::file_path_as_absolute(x)
  }))


  # additional check for STARsolo
  check_STARsolo_output <- function(x) {
    files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
    test <- list.files(x) %in% c(files, paste0(files, ".gz"))
    length(test) != 0 && all(test)
  }


  if (!"WellBC" %in% colnames(sample_file)) {
    stop(paste0("STARsolo output is used. The sample annotation must ",
                "contain 'WellBC' column providing cell barcodes."))
  }

  mapping_dir <- unique(mapping_dir)
  all_STARsolo_files_present <- all(
    unlist(
      lapply(mapping_dir, function(x) {
        check_STARsolo_output(x)
      })
    )
  )
  if (!all_STARsolo_files_present) {
    stop(paste0("Could not find files: 'barcodes', 'features' and 'matrix'",
                " for STARsolo output. Please check 'mappingDir' argument."))
  }

  message("- 'SampleFileName' column - OK")



  list(sample_expression_files = mapping_dir)
}

create_exprs_matrix <- function(exprs_matrix_path, exprs_file_paths,
                                output, measure, col_names, cell_barcodes) {

  read_matrix <- Read10X(data_dir = exprs_file_paths, gene_column = 1)
  # keep index of feature names containing "_" because Seurat
  #changes them to "-" and they no longer match with fdata[, "gene_id"]
  idx <- grep("_", rownames(read_matrix))

  requireNamespace("Seurat")
  seurat_object <- Seurat::CreateSeuratObject(counts = read_matrix)

  exprs_matrix <- as.matrix(seurat_object[['RNA']]$counts)
  # replace "-" with "_" for features with "_" 
  # before converting to Seurat object
  rownames(exprs_matrix)[idx] <- gsub("-", "_", rownames(exprs_matrix)[idx])
  requireNamespace("stringr")
  exprs_matrix <- exprs_matrix[, stringr::str_detect(colnames(exprs_matrix),
                                  paste(cell_barcodes, collapse = "|"))]


  # check if rownames are ENSEMBL and remove version suffix
  isENSEMBL <- all(grepl("ENS", rownames(exprs_matrix)))
  if (isENSEMBL) {
    # do not use gsub("(.+)[.]\\d+", "\\1", rownames(exprs_matrix)),
    # so that ENS000000.1_PAR_Y can be kept
    rownames(exprs_matrix) <- gsub("\\.\\d+$", "", rownames(exprs_matrix))
  }


  colnames(exprs_matrix) <- col_names

  exprs_matrix
}

create_eset <- function(feature_annotation_path,
                        sample_annotation_path,
                        mapping_dir,
                        barcodes,
                        output_path,
                        pool_name,
                        exprs_matrix_path = NULL,
                        path = NULL,
                        add_eset_annotation = NULL) {
  if (!file.exists(feature_annotation_path)) {
    stop("Could not find feature annotation at '", feature_annotation_path, "'")
  }

  if (!file.exists(sample_annotation_path)) {
    stop("Could not find sample annotation at '", sample_annotation_path, "'")
  }

  if(!is.null(exprs_matrix_path)) {
    if(!file.exists(exprs_matrix_path)) {
      stop("Could not find expression matrix at '", exprs_matrix_path, "'")
    }
  }

  if(!is.null(path)) {
    if(!dir.exists(path)) {
      stop("Provided 'path': '", path, "' does not exist.")
    }
  }

  ##### Import annotation files #####
  message("Importing feature annotation")
  fdata_file <- read.table(feature_annotation_path, header = TRUE,
                           sep = "\t", quote = "\"",
                           comment.char = "", stringsAsFactors = FALSE)

  # for backwards compatibility
  if("ENSEMBL" %in% colnames(fdata_file) && !all(grepl("ENS", fdata_file[, "ENSEMBL"])) & !"gene_id" %in% colnames(fdata_file)) {
    colnames(fdata_file)[which(colnames(fdata_file) == "ENSEMBL")] <- "gene_id"
  }

  # Check gene annotation
  if(!"gene_id" %in% colnames(fdata_file))
    stop("'gene_id' column with unique feature identifiers must be present in 'feature_annotation_path'.")

  # check if duplicated ids are present
  if(any(duplicated(fdata_file$gene_id)))
    stop("Duplicated features ids are not allowed. Please check the 'gene_id' column in 'feature_annotation_path'.")

  message("Importing sample annotation")
  sample_file <- read.table(sample_annotation_path, header = TRUE,
                            sep = "\t", quote = "\"",
                            comment.char = "", stringsAsFactors = FALSE)
  # Check sample annotation
  check_sample_file_list <- check_sample_file(mapping_dir = mapping_dir,
                                              sample_file = sample_file)
  output <- "STARsolo"
  measure <- "counts"
  sample_expression_files <- check_sample_file_list$sample_expression_files

  ##### Create phenodata #####
  pdata_eset <- create_pdata(sample_file = sample_file, pool_name = pool_name,
                             barcodes = barcodes)

  ##### Create expression matrix #####
  message("Creating expression matrix")

  exprs_matrix_eset <- create_exprs_matrix(
    exprs_matrix_path = exprs_matrix_path,
    exprs_file_paths = sample_expression_files,
    output = output,
    measure = measure,
    col_names = rownames(pdata_eset),
    cell_barcodes = barcodes
  )


  ##### Create featuredata #####
  message("Creating feature data")

  fdata_eset <- fdata_file
  rownames(fdata_eset) <- fdata_eset[, "gene_id"]

  # intersect features between exprs matrix and fdata
  feature_files <- match_features(exprs_matrix = exprs_matrix_eset,
                                  fdata = fdata_eset)

  fdata_eset <- feature_files$fdata
  exprs_matrix_eset <- feature_files$exprs_matrix

  ##### Create eSet #####
  message("Creating eset")

  if (nrow(pdata_eset) != ncol(exprs_matrix_eset)) {
    stop("nrow(pData) and ncol(exprsMatrix) differ")
  }

  if (nrow(fdata_eset) != nrow(exprs_matrix_eset)) {
    stop("nrow(fData) and nrow(exprsMatrix) differ")
  }

  if (!all(rownames(pdata_eset) == colnames(exprs_matrix_eset))) {
    stop("rownames(pData) and colnames(exprsMatrix) differ")
  }

  if (!all(rownames(fdata_eset) == rownames(exprs_matrix_eset))) {
    stop("rownames(fData) and rownames(exprsMatrix) differ")
  }

  if (!inherits(exprs_matrix_eset, "matrix")) {
    stop("exprsMatrix must be of class 'matrix'")
  }



  additional_info <- paste0("Additional information about eSet \n",
                            "  Expression matrix created from ",
                            output, " output. \n",
                            "  Expression matrix contains non-transformed ",
                            ifelse(output %in% c("STAR", "STARsolo"),
                                   "counts",
                                   ifelse(measure == "expected_count",
                                          "counts", measure)), ".")


  if (isTRUE(!is.null(add_eset_annotation) &
               is.character(add_eset_annotation))) {
    additional_info <- paste0(additional_info, "\n", "  ", add_eset_annotation)
  }

  fdata_eset <- new("AnnotatedDataFrame", data = fdata_eset)
  pdata_eset <- new("AnnotatedDataFrame", data = pdata_eset)

  requireNamespace("Biobase")
  eset <- Biobase::ExpressionSet(assayData = exprs_matrix_eset,
                                  phenoData = pdata_eset,
                                  featureData = fdata_eset,
                                  annotation = additional_info)

  eset <- eset[, colSums(exprs(eset)) != 0]
  saveRDS(eset, file = output_path)

  message(paste0("eset created succesfully for ", ncol(eset),
                 " samples and ", nrow(eset),
                 " genes and saved at ", output_path, ".")) 

  eset
}


p_data_file <- par$pDataFile
f_data_file <- par$fDataFile
pool_name <- par$poolName
mapping_dir <- lapply(par$mappingDir,
                      \(x) file.path(x, "Solo.out", "Gene", "raw"))

get_barcode_from_mapping_dir <- function(raw_dir) {
  barcodes_file <- file.path(raw_dir, "barcodes.tsv")
  if (!file.exists(barcodes_file)) {
    stop(paste0("Expected the 'Solo.out/Gene/raw' directory at ",
                raw_dir, " to contain a 'barcodes.tsv' file."))
  }
  barcodes <- readLines(barcodes_file)
  if (length(barcodes) != 1) {
    stop(paste0("A single STAR Solo folder should only have ",
                "mapped one (1) barcode, but found '",
                length(barcodes), "'for mapping directory ", raw_dir))
  }
  return(barcodes)
}

barcodes <- lapply(mapping_dir, get_barcode_from_mapping_dir)

print(paste0("mappingDir: ", mapping_dir))
print(paste0("pDataFile: ", p_data_file))
print(paste0("fDataFile: ", f_data_file))
print(paste0("poolName: ", pool_name))
print(paste0("barcodes: ", barcodes))



# CREATE ESET WITH RAW UMI COUNTS

eset <- create_eset(feature_annotation_path = f_data_file,
                    sample_annotation_path = p_data_file,
                    mapping_dir = mapping_dir,
                    barcodes = barcodes,
                    output_path = par$output,
                    pool_name = pool_name,
                    path = NULL,
                    exprs_matrix_path = NULL)