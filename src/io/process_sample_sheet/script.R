library(data.table)

## VIASH START
par <- list(
  "id" = "run_id",
  "sample_sheet" = "sample_sheet1.csv",
  # "sample_sheet" = "sample_sheet2.csv",
  "project_name" = "project_name",
  "instrument" = "Foo",
  "instrument_type" = "Bar",
  "run_info" = "run_info.tsv",
  "plate_info" = "plate_info.tsv"
)
## VIASH END

md_header <- function(x) {
  cat(paste0("\n", x, "\n\n"))
}

md_item <- function(x) {
  cat(x, "\n")
}

md_block <- function(x) {
  print(paste0(x, "\n"))
}

md_header("# Process sample sheet")

id <- par$id
sample_sheet_file <- par$sample_sheet
project <- par$project_name
instrument <- par$instrument
type <- par$instrument_type
project_name <- par$project_name
run_info_file <- par$run_info
plate_info_file <- par$plate_info

md_header(paste("## Parsing sample sheet file:", sample_sheet_file))

sample_sheet_lines <- readLines(sample_sheet_file)

# Where does the sample section start?
start_data_line <- which(grepl("^\\[[dD]ata\\]|^\\[Cloud_Data\\]", sample_sheet_lines))
if (length(start_data_line) == 0) {
  stop("No data section found in the sampleSheet")
}
md_item(paste("- Sample information starts at: ", start_data_line))

# Read the sample sheet csv section
plate_data <- data.table::fread(
  sample_sheet_file,
  skip = start_data_line,
  colClasses = "character",
  header = TRUE
)

md_item("- Sample sheet data:")
print(plate_data)

md_header("## Check if project name is provided in the sample sheet")
if (!"ProjectName" %in% names(plate_data)) {
  if (is.null(project_name)) {
    stop("ProjectName column is not provided in sampleSheet or parameter")
  } else {
    md_header("### Setting project name from input parameter")
    plate_data$ProjectName <- rep(project_name, length(plate_data$Sample_ID))
  }
}

# TODO: Should this check be here?!
# plate_data <- plate_data[which(ProjectName == project)]

md_header("## Check if instrument and type are provided in the sample sheet")
cat("If not, the user can specify them as parameters", "\n\n")

instrument_info <- plate_data$Instrument
type_info <- plate_data$Type

md_header("### Instrument and type from sample sheet")

md_item(paste("- Instrument info from samplesheet: ", instrument_info))
md_item(paste("- Instrument type info from samplesheet: ", type_info))

if (length(unique(instrument_info)) > 1 || length(unique(type_info)) > 1) {
  stop("Instrument and/or type are not unique in the sampleSheet")
} else {
  instrument_info <- unique(instrument_info)
  type_info <- unique(type_info)
}

if (length(instrument_info) != 1 || length(type_info) != 1) {
  if (is.null(instrument) || is.null(type)) {
    stop("Instrument and/or type are not specified in sampleSheet nor parameters")
  } else {
    md_header("### Setting instrument and type from input parameters")
    instrument_info <- instrument
    type_info <- type
  }
}

md_item(paste("- Instrument info: ", instrument_info))
md_item(paste("- Instrument type info: ", type_info))

md_header("## Write run metadata to file and add instrument and type information")
run_data <- data.table(
  id = id,
  instrument = instrument_info,
  type = type_info
)

print(run_data)

data.table::fwrite(
  run_data,
  file = run_info_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

md_header("## Validate sample metadata and write to file")

if (!all(c("Sample_ID", "I7_Index_ID") %in% names(plate_data))) {
  stop("Sample_ID and I7_Index_ID columns are not present in the sampleSheet")
} else {
  setnames(plate_data, old = c("Sample_ID", "I7_Index_ID"), new = c("PoolName", "IDT"))

  print(plate_data)

  if ("index" %in% names(plate_data)) {
    data.table::fwrite(
      plate_data[, .(PoolName, ProjectName, IDT, index)],
      file = plate_info_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  } else {
    data.table::fwrite(
      plate_data[, .(PoolName, ProjectName, IDT)],
      file = plate_info_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }
}
