library(data.table)

## VIASH START
par <- list(
  "id" = "run_id",
  # "sample_sheet" = "sample_sheet1.csv",
  "sample_sheet" = "sample_sheet2.csv",
  "project_name" = "project_name",
  # "instrument" = "Foo",
  # "instrument_type" = "Bar",
  "run_info" = "run_info.tsv",
  "plate_info" = "plate_info.tsv"
)
## VIASH END

id <- par$id
sample_sheet_file <- par$sample_sheet
project <- par$project_name
instrument <- par$instrument
type <- par$instrument_type
run_info_file <- par$run_info
plate_info_file <- par$plate_info

cat(">> Parsing sample sheet file: ", sample_sheet_file, "\n")

sample_sheet_lines <- readLines(sample_sheet_file)

# Where does the sample section start?
start_data_line <- which(grepl("^\\[[dD]ata\\]|^\\[Cloud_Data\\]", sample_sheet_lines))
if (length(start_data_line) == 0) {
  stop("No data section found in the sampleSheet")
}
cat(paste(">>> Sample information starts at: ", start_data_line), "\n")

# Read the sample sheet csv section
plate_data <- data.table::fread(
  sample_sheet_file,
  skip = start_data_line,
  colClasses = "character",
  header = TRUE
)

cat(">>> Sample sheet data:", "\n")
print(plate_data)

# TODO: Should this check be here?!
# plate_data <- plate_data[which(project_name == project)]

cat(">> Check if instrument and type are provided in the sample sheet", "\n")
cat("   If not, the user can specify them as parameters", "\n")
instrument_info <- plate_data$Instrument
type_info <- plate_data$Type
cat(paste(">>> Instrument info from samplesheet: ", instrument_info), "\n")
cat(paste(">>> Instrument type info from samplesheet: ", type_info), "\n")

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
    cat(">>> Setting instrument and type from input parameters", "\n")
    instrument_info <- instrument
    type_info <- type
  }
}

cat(paste(">>> Instrument info: ", instrument_info), "\n")
cat(paste(">>> Instrument type info: ", type_info), "\n")

cat(">> Write run metadata to file and add instrument and type information", "\n")
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

cat(">> Validate sample metadata and write to file", "\n")

if (!all(c("Sample_ID", "I7_Index_ID") %in% names(plate_data))) {
  stop("Sample_ID and I7_Index_ID columns are not present in the sampleSheet")
}

setnames(plate_data, old = c("Sample_ID", "I7_Index_ID"), new = c("PoolName", "IDT"))

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
