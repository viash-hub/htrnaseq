#!/bin/bash

## VIASH START
# par_input="Homo_sapiens/L001/,Homo_sapiens/L002/,Homo_sapiens/L003/,Homo_sapiens/L004/"
# par_barcodes="AACAAGGTAC,GTCTCGAGTG,GTTGAATTGG"
# par_output="output-test"
par_input="E090713_TOX15991_S1__001.cutadapt.outputDir"
par_poolName="test_pool"
par_barcodes="AACAAGGTAC,GTCTCGAGTG,GTTGAATTGG"
par_output="output-test"
par_genomeDir="star"
par_wellBarcodesLength=10
par_umiLength=10
par_limitBAMsortRAM="10000000000"
meta_cpus=2
## VIASH END

# set -ex
if [ -z $par_output ]; then
  par_output=.
else
  mkdir -p "$par_output"
fi

barcodes=$(echo $par_barcodes | tr "," " ")

# Load reference genome
STAR --genomeLoad LoadAndExit --genomeDir $par_genomeDir

function _run() {

  local barcode="$1"
  local par_wellBarcodeLength="$2"
  local par_UMIlength="$3"
  local par_input="$4"
  local par_poolName="$5"
  local par_output="$6"
  local par_genomeDir="$7"
  local par_limitBAMsortRAM="$8"
  local par_runThreadN="$9"

  local par_UMIstart=$(($par_wellBarcodeLength + 1))

  # Prepare the input for STAR
  # This needs to be of the format:
  #   'R2 R1'
  local R1=""
  local R2=""

  for lane in $lanes; do
    R1="$lane/${barcode}_R1_001.fastq"
    R2="$lane/${barcode}_R2_001.fastq"
    R1s="$R1s $R1"
    R2s="$R2s $R2"
  done

  local input_R2=`echo $R2s | tr ' ' ','`
  local input_R1=`echo $R1s | tr ' ' ','`
  local star_readFilesIn="$input_R2 $input_R1"

  echo
  echo "Processing $barcode for $par_poolName"
  echo "For the following inputs (lanes):"
  echo "$par_input" | tr ',' '\n'
  echo
  echo "$star_readFilesIn"

  # Write <barcode> to <barcode>.txt
  echo "$barcode" > "$barcode.txt"

  local dir="$par_output/$par_poolName/$barcode/"
  mkdir -p "$dir"

  STAR \
    --readFilesIn $star_readFilesIn \
    --soloType Droplet \
    --quantMode GeneCounts \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM $par_limitBAMsortRAM \
    --runThreadN $par_runThreadN \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --soloCBstart 1 \
    --soloCBlen $par_wellBarcodeLength \
    --soloUMIstart $par_UMIstart \
    --soloUMIlen $par_UMIlength \
    --soloBarcodeReadLength 0 \
    --soloStrand Unstranded \
    --soloFeatures Gene \
    --genomeDir "$par_genomeDir" \
    --outReadsUnmapped Fastx \
    --outSAMunmapped Within \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
    --soloCBwhitelist "$barcode.txt" \
    --outFileNamePrefix $dir

}

# Export the function - requires bash
export -f _run

# Run the jobs in parallel
# We ask GNU parallel to retry failed jobs. This might have an impact on the derivation of the number of errors later
export PARALLEL_SHELL="/bin/bash"
parallel --jobs 80% --memfree 2G --retry-failed --retries 4 --halt now,fail=1 --joblog execution_log.txt '_run {}' ::: $barcodes ::: $par_wellBarcodesLength ::: $par_umiLength ::: $par_input ::: $par_poolName ::: $par_output ::: $par_genomeDir ::: $par_limitBAMsortRAM ::: $par_runThreadN

# Unload reference
STAR --genomeLoad Remove --genomeDir $par_genomeDir

# Handle errors
errors=$(cat execution_log.txt | grep -v Exitval | cut -f7 | grep -vw 0 | wc -l)
if ((errors>0)); then
  echo ""
  echo "=================================================================="
  echo ""
  echo "!!! $errors error(s) occurred"
  echo ""
  cat execution_log.txt
  echo ""
  echo "=================================================================="
  echo ""
  exit 1
else
  echo ""
  echo "=================================================================="
  echo ""
  echo "Mapping went fine, zero errors occurred"
  echo ""
  echo "=================================================================="
  echo ""
fi






