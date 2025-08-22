#!/bin/bash

set -eo pipefail

echo "Publishing results to multiple output directories"

# Create output directories
echo "Creating output directories..."
mkdir -p "$par_star_output_dir" && echo "$par_star_output_dir created"
mkdir -p "$par_qc_metrics_dir" && echo "$par_qc_metrics_dir created"  
mkdir -p "$par_eset_dir" && echo "$par_eset_dir created"
mkdir -p "$par_data_dir" && echo "$par_data_dir created"
mkdir -p "$par_reports_dir" && echo "$par_reports_dir created"

echo
echo "Copying STAR output files..."
IFS=";" read -ra star_output <<<$par_star_output
for i in "${star_output[@]}"; do
  echo "Copying $i to $par_star_output_dir/"
  cp -rL "$i" "$par_star_output_dir/"
done

echo
echo "Copying QC metrics files..."
IFS=";" read -ra nrReadsNrGenesPerChrom <<<$par_nrReadsNrGenesPerChrom
for i in "${nrReadsNrGenesPerChrom[@]}"; do
  echo "Copying $i to $par_qc_metrics_dir/"
  cp -rL "$i" "$par_qc_metrics_dir/"
done

IFS=";" read -ra star_qc_metrics <<<$par_star_qc_metrics
for i in "${star_qc_metrics[@]}"; do
  echo "Copying $i to $par_qc_metrics_dir/"
  cp -rL "$i" "$par_qc_metrics_dir/"
done

echo
echo "Copying eset files..."
IFS=";" read -ra eset <<<$par_eset
for i in "${eset[@]}"; do
  echo "Copying $i to $par_eset_dir/"
  cp -rL "$i" "$par_eset_dir/"
done

echo
echo "Copying data files..."
IFS=";" read -ra f_data <<<$par_f_data
for i in "${f_data[@]}"; do
  echo "Copying $i to $par_data_dir/"
  cp -rL "$i" "$par_data_dir/"
done

IFS=";" read -ra p_data <<<$par_p_data
for i in "${p_data[@]}"; do
  echo "Copying $i to $par_data_dir/"
  cp -rL "$i" "$par_data_dir/"
done

echo
echo "Copying reports and run parameters..."
echo "Copying $par_html_report to $par_reports_dir/"
cp -rL "$par_html_report" "$par_reports_dir/"

echo "Copying $par_run_params to $par_reports_dir/"
cp -rL "$par_run_params" "$par_reports_dir/"

echo
echo "Publishing completed successfully!"
