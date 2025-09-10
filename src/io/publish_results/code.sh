#!/bin/bash

set -eo pipefail

echo "Publishing results to multiple output directories"

# Create output directories for multiple files
echo "Creating output directories..."

path_pars=(
  par_star_output_dir
  par_nrReadsNrGenesPerChrom_dir
  par_star_qc_metrics_dir
  par_eset_dir
  par_f_data_dir
  par_p_data_dir
  par_html_report_output
  par_run_params_output
)

for par in ${path_pars[@]}; do
    curr_val="${!par}"
    new_value=$(realpath --canonicalize-missing "$curr_val")
    declare -g "$par=$new_value"
done

mkdir -p "$par_star_output_dir" && echo "$par_star_output_dir created"
mkdir -p "$par_nrReadsNrGenesPerChrom_dir" && echo "$par_nrReadsNrGenesPerChrom_dir created"  
mkdir -p "$par_star_qc_metrics_dir" && echo "$par_star_qc_metrics_dir created"
mkdir -p "$par_eset_dir" && echo "$par_eset_dir created"
mkdir -p "$par_f_data_dir" && echo "$par_f_data_dir created"
mkdir -p "$par_p_data_dir" && echo "$par_p_data_dir created"

echo
echo "Copying STAR output files..."
IFS=";" read -ra star_output <<<$par_star_output
for i in "${star_output[@]}"; do
  echo "Copying $i to $par_star_output_dir/"
  cp -a --keep-directory-symlink "$i" "$par_star_output_dir/"
done

echo
echo "Copying nrReadsNrGenesPerChrom files..."
IFS=";" read -ra nrReadsNrGenesPerChrom <<<$par_nrReadsNrGenesPerChrom
for i in "${nrReadsNrGenesPerChrom[@]}"; do
  echo "Copying $i to $par_nrReadsNrGenesPerChrom_dir/"
  cp -a --keep-directory-symlink "$i" "$par_nrReadsNrGenesPerChrom_dir/"
done

echo
echo "Copying STAR QC metrics files..."
IFS=";" read -ra star_qc_metrics <<<$par_star_qc_metrics
for i in "${star_qc_metrics[@]}"; do
  echo "Copying $i to $par_star_qc_metrics_dir/"
  cp -a --keep-directory-symlink "$i" "$par_star_qc_metrics_dir/"
done

echo
echo "Copying eset files..."
IFS=";" read -ra eset <<<$par_eset
for i in "${eset[@]}"; do
  echo "Copying $i to $par_eset_dir/"
  cp -a --keep-directory-symlink "$i" "$par_eset_dir/"
done

echo
echo "Copying f_data files..."
IFS=";" read -ra f_data <<<$par_f_data
for i in "${f_data[@]}"; do
  echo "Copying $i to $par_f_data_dir/"
  cp -a --keep-directory-symlink "$i" "$par_f_data_dir/"
done

echo
echo "Copying p_data files..."
IFS=";" read -ra p_data <<<$par_p_data
for i in "${p_data[@]}"; do
  echo "Copying $i to $par_p_data_dir/"
  cp -a --keep-directory-symlink "$i" "$par_p_data_dir/"
done

echo
echo "Copying single files directly..."
mkdir -p $(dirname "$par_html_report_output")
echo "Copying $par_html_report to $par_html_report_output"
cp -a --keep-directory-symlink "$par_html_report" "$par_html_report_output"

echo "Copying $par_run_params to $par_run_params_output"
mkdir -p $(dirname "$par_run_params_output")
cp -a --keep-directory-symlink "$par_run_params" "$par_run_params_output"

echo
echo "Publishing completed successfully!"
