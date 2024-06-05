library(data.table)

## VIASH START
par <- list(
  "id" = "id",
  "sample_sheet" = "sample_sheet1.csv",
  # "sample_sheet" = "sample_sheet2.csv",
  "project_name" = "project_name",
  "instrument" = "Foo",
  "instrument_type" = "Bar"
)
## VIASH END

id <- par$id
sample_sheet_file <- par$sample_sheet
project <- par$project_name
instrument <- par$instrument
type <- par$instrument_type
# run_info_file <- par$run_info
# plate_info_file <- par$plate_info

cat(">> Parsing sample sheet file: ", sample_sheet_file)

sample_sheet_lines <- readLines(sample_sheet_file)

# Where does the sample section start?
start_data_line <- which(grepl("^\\[Data\\]", sample_sheet_lines))
cat(paste(">>> Sample information starts at: ", start_data_line))

# Read the sample sheet csv section
plate_data <- data.table::fread(
  sample_sheet_file,
  skip = start_data_line,
  colClasses = "character",
  header = TRUE
)

print(">>> Sample sheet data:")
print(plate_data)

# TODO: Should this check be here?!
# plate_data <- plate_data[which(project_name == project)]

print(">> Check if instrument and type are provided in the sample sheet")
print(">> If not, the user can specify them as parameters as well")
instrument_info <- plate_data$Instrument
type_info <- plate_data$type

if (length(instrument_info) > 1 || length(type_info) > 1) {
  stop("Instrument and/or type are not uniquely specified in the sampleSheet")
} else {
  instrument_info <- unique(instrument_info)
  type_info <- unique(type_info)
}

if (length(instrument_info) != 1 || length(type_info) != 1) {
  if (is.null(instrument) || is.null(type)) {
    stop("Instrument and/or type are not specified in the sampleSheet nor as input parameters")
  } else {
    instrument_info <- instrument
    type_info <- type
  }
}

cat(paste(">>> Instrument info: ", instrument_info))
cat(paste(">>> Instrument type info: ", type_info))

cat(">> Write the run and plate info to files and add instrument and type information")
runData <- data.table(
  id = id,
  instrument = instrument_info,
  type = type_info
)

print(runData)



#
#
# setnames(plateData, old = c("Sample_ID", "I7_Index_ID"), new = c("PoolName", "IDT"))
#
# if ("index" %in% names(plateData)) {
#   data.table::fwrite(plateData[, .(PoolName, ProjectName, IDT, index)], file = plateInfoFile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# } else {
#   data.table::fwrite(plateData[, .(PoolName, ProjectName, IDT)], file = plateInfoFile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# }
# data.table::fwrite(runData, file = runInfoFile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
