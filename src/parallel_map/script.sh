#!/bin/bash

## VIASH START
par_input_r1="work/2c/5b8b3a2dd4a988b8838e3f72d38a37/_viash_par/input_r1_1/two__ACACCGAATT.concat_text_r1.output.txt"
par_input_r2="work/2c/5b8b3a2dd4a988b8838e3f72d38a37/_viash_par/input_r2_1/two__ACACCGAATT.concat_text_r2.output.txt"
par_barcodes="ACACCGAATT;GGCTATTGAT"
par_output="./*"
par_genomeDir="star"
par_wellBarcodesLength=10
par_umiLength=10
par_limitBAMsortRAM="10000000000"
meta_cpus=2
par_runThreadN=1
## VIASH END

set -eo pipefail

# Check if wildcard character is present in output folder template
printf "Checking if output folder template ($par_output) contains a single wildcard character '*'. "
output_glob_character="${par_output//[^\*]}"
if [[ "${#output_glob_character}" -ne "1" ]]; then
  echo "The value for --output must contain exactly one '*' character. Exiting..."
  exit 1
else
  echo "Done, wildcard character found!"
fi

# Split the delimited strings into arrays
IFS=';' read -r -a barcodes <<< "$par_barcodes"
IFS=';' read -r -a input_r1 <<< "$par_input_r1"
IFS=';' read -r -a input_r2 <<< "$par_input_r2"

# Check that the number of values provided for the barcodes and the fastq files are the same.
num_barcodes="${#barcodes[@]}"
num_r1_inputs="${#input_r1[@]}"
num_r2_inputs="${#input_r2[@]}"

if [ ! "$num_barcodes" -eq "$num_r1_inputs" ] || [ ! "$num_r1_inputs" -eq "$num_r2_inputs" ]; then
  echo "The number of values for arguments 'barcodes' ($num_barcodes), "\
        "'input_r1' ($num_r1_inputs) and 'input_r2' ($num_r2_inputs) "\
        "should be the same, and their order should match."
  exit 1
else
  echo "Checked if length of barcodes input ($num_barcodes) is "\
       "the same as R1 reads ($num_r1_inputs) and R2 reads "\
       "($num_r2_inputs). Seems OK!"
fi

# Function to test for unique values in array
function arrayContainsUniqueValues {
  # Pass the argument by reference
  local -n arr=$1
  # Create a temporary associative array
  # in order to use its uniqueness of keys
  # 'declare' in a function is automatically local
  declare -A uniq_tmp
  for item in "${arr[@]}"; do
    uniq_tmp[$item]=0 # assigning a placeholder
  done
  local unique_array_values=(${!uniq_tmp[@]})
  if [ "${#unique_array_values[@]}" -eq "${#arr[@]}" ]; then
    return
  fi
  false
}
arrayContainsUniqueValues barcodes
is_array_unique_exit_code=$?
if ! (exit $is_array_unique_exit_code); then 
  echo "The provided barcodes should be unique!"
  echo "Values: $par_barcodes"
  exit 1
fi

# Define the function that will be used to run a single job
function _run() {
  local par_wellBarcodeLength="$1"
  local par_UMIlength="$2"
  local par_output="$3"
  local par_genomeDir="$4"
  local par_limitBAMsortRAM="$5"
  local par_runThreadN="$6"
  local barcode="$7"
  local input_R1="$8"
  local input_R2="$9"
  local par_UMIstart=$(($par_wellBarcodeLength + 1))

  set -eo pipefail

  echo <<-EOF
    Processing $barcode
    For the following inputs (lanes):
    "$star_readFilesIn
	EOF

  echo "Writing barcode '$barcode' to $barcode.txt and using it as input".
  # Note that there is no possible conflict between jobs here
  # because the barcodes are unique (and the barcode is part of the name
  # of the file).
  echo "$barcode" > "$barcode.txt"

  local dir="${par_output//\*/$barcode}/"
  echo "Setting output for barcode '$barcode' to '$dir'."
  mkdir -p "$dir"

  # check if files are compressed
  local TMPDIR=$(mktemp -d "$meta_temp_dir/parallel_map-$barcode-XXXXXX")
  function clean_up {
    [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
  }
  trap clean_up RETURN

  # Decompress the input files when needed
  # NOTE: for some reason, using STAR's --readFilesCommand does not always work
  # This might be because STAR creates fifo files (see https://man7.org/linux/man-pages/man7/fifo.7.html)
  # and this requires a filesystem that supports this. Another cause might be that the input files
  # are symlinks. When testing this, using '--readFilesCommand "zcat"' 
  # always produced empty BAM files, but also a succesfull exit code (0) so the problem is not reported.
  # However, the logs showed the following error: "gzip -: unexpected end of file".

  function is_gzipped {
    printf "Checking if input '$1' (barcode '$barcode') is gzipped... "
    if file "$1" | grep -q 'gzip'; then
      echo "Done, detected compressed file."
      return
    fi
    echo "Done, file does not need decompression."
    false
  }
  
  # Resolve symbolic links to actual file paths
  input_R1=$(realpath $input_R1)
  input_R2=$(realpath $input_R2)

  if is_gzipped $input_R1; then
    local compressed_file_name_r1="$(basename -- $input_R1)"
    local uncompressed_file_r1="$TMPDIR/${compressed_file_name_r1%.gz}"
    printf "Unpacking input to $uncompressed_file_r1... "
    zcat "$input_R1" > "$uncompressed_file_r1"
    echo "Decompression done."
  else
    local uncompressed_file_r1="$input_R1"
  fi

  if is_gzipped $input_R2; then
    local compressed_file_name_r2="$(basename -- $input_R2)"
    local uncompressed_file_r2="$TMPDIR/${compressed_file_name_r2%.gz}"
    printf "Unpacking input to $uncompressed_file_r2... "
    zcat "$input_R2" > "$uncompressed_file_r2"
    echo "Decompression done."
  else
    local uncompressed_file_r2="$input_R2"
  fi

  local n_input_lines_r1=$(wc -l < "$uncompressed_file_r1")
  local n_input_lines_r2=$(wc -l < "$uncompressed_file_r2")

  printf "Checking if length of input file mates match. "
  if (( $n_input_lines_r1 != n_input_lines_r2 )); then
    echo "The length of file $input_R1 ($n_input_lines_r1) does not match with $input_R2 ($n_input_lines_r2)"
    return 1
  else
    echo "Seems OK, $n_input_lines_r1 input lines."
  fi
  echo "Starting STAR for barcode '$barcode'"
  # soloType 'Droplet' is the same as 'CB_UMI_Simple': one UMI and one cell barcode of fixed length. 
  # By default in this mode, STAR will look for the cell barcode and the UMI int the last files specified with --readFilesIn
  # So we need to specify R2 first and R1 second, because R1 contains the barcode and UMI.
  # Also, you might be tempted to use '--soloBarcodeMate 1' to alter this behavior, but this requires the clipping
  # the barcode from this mate by specifying --clip5pNbases and/or --clip3pNbases, which we do not want to do.
  STAR \
    --readFilesIn "$uncompressed_file_r2" "$uncompressed_file_r1" \
    --soloType Droplet \
    --quantMode GeneCounts \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM "$par_limitBAMsortRAM" \
    --runThreadN "$par_runThreadN" \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --soloCBstart 1 \
    --readFilesType "Fastx" \
    --soloCBlen "$par_wellBarcodeLength" \
    --soloUMIstart "$par_UMIstart" \
    --soloUMIlen "$par_UMIlength" \
    --soloBarcodeReadLength 0 \
    --soloStrand Unstranded \
    --soloFeatures Gene \
    --genomeDir "$par_genomeDir" \
    --outReadsUnmapped Fastx \
    --outSAMunmapped Within \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
    --soloCBwhitelist "$barcode.txt" \
    --outFileNamePrefix "$dir" \
    --outTmpDir "$TMPDIR/STARtemp/"

  printf "Done running STAR. "
  # Check if the number of processed reads is equal to the number of input reads
  local n_input_reads=$(($n_input_lines_r1 / 4))
  local nr_output_reads=$(grep -Po "Number\ of\ input\ reads \\|\W*\K\d+" "$dir/Log.final.out")
  if (( $nr_output_reads != $n_input_reads )); then
    echo "Not all input reads were processed for barcode $barcode."
    return 1
  else
    echo "Processed $nr_output_reads reads for barcode $barcode".
  fi

  printf "Making sure that the output has the proper permissions."
  find "$dir" -type d -exec chmod o+x {} \;
  chmod -R o+r "$dir"
  echo "Done"
}

# Export the function - requires bash
export -f _run

# Load reference genome
echo "Loading reference genome"
STAR --genomeLoad LoadAndExit --genomeDir "$par_genomeDir"

# Run the concurrent jobs using GNU parallel

# Make sure that parallel uses the correct shell
export PARALLEL_SHELL="/bin/bash"

# Some notes:
#   --halt now,fail=1: instruct parallel to exit when a job has failed and kill remaining running jobs.
#   
# ::: is a special syntax for GNU parallel to delineate inputs
# If multiple ::: are given, each group will be treated as an input source, and all combinations of input
# sources will be generated. E.g. ::: 1 2 ::: a b c will result in the combinations (1,a) (1,b) (1,c) (2,a) (2,b) (2,c)
# The delimiter :::+ (note the extra '+') links the argument to the previous argument, and one argument from each of the input
# sources will be read.
parallel_cmd=("parallel" "--jobs" "80%" "--verbose" "--memfree" "2G"
              "--tmpdir" "$meta_temp_dir"
              "--retry-failed" "--retries" "4" "--halt" "soon,fail=1"
              "--joblog" "$par_joblog" "_run" "{}")

# Arguments for which there is one value, so these will not create extra jobs
parallel_cmd+=(":::" "$par_wellBarcodesLength" ":::" "$par_umiLength" ":::" "$par_output" ":::" "$par_genomeDir" ":::" "$par_limitBAMsortRAM" ":::" "$par_runThreadN")

# Argument which in fact will cause extra jobs to be spawned, per job one item from each argument will be selected
# Thus, these argument lists should have the same length.
parallel_cmd+=(":::" "${barcodes[@]}" ":::+" "${input_r1[@]}" ":::+" "${input_r2[@]}")

set +eo pipefail
"${parallel_cmd[@]}"
exit_code=$?
set -eo pipefail

echo "GNU parallel finished!"

# Unload reference
printf "Unloading reference genome. "
STAR --genomeLoad Remove --genomeDir "$par_genomeDir"
echo "Done!"

# Exit code from GNU parallel:
# If fail=1 is used, the exit status will be the exit status of the failing job.
echo "Checking exit code"
if ((exit_code>0)); then
  # Note that the ending HERE must be indented with TAB characters (not spaces)
  # in order to remove leading indentation
  MESSAGE=$(
    cat <<-HERE
    ==================================================================

    !!! An error occurred for one of the jobs.
    Exit code of the failing job: $exit_code.

    %s

    ==================================================================

		HERE
  )
  printf "$MESSAGE" "$(<$par_joblog)"
  exit 1
else
  cat <<-HERE
  ==================================================================

  Mapping went fine (exit code '$exit_code'), zero errors occurred

  ==================================================================
	HERE

fi






