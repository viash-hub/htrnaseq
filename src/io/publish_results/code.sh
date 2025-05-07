#!/bin/bash

echo "Publishing $par_input -> $par_output"

echo
echo "Creating directory if it does not exist:"
mkdir -p "$par_output" && echo "$par_output created"

echo
echo "Copying files..."
IFS=";" read -ra star_output <<<$par_star_output
IFS=";" read -ra nrReadsNrGenesPerChrom <<<$par_nrReadsNrGenesPerChrom
IFS=";" read -ra star_qc_metrics <<<$par_star_qc_metrics
IFS=";" read -ra eset <<<$par_eset
IFS=";" read -ra f_data <<<$par_f_data
IFS=";" read -ra p_data <<<$par_p_data

for i in "${star_output[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${nrReadsNrGenesPerChrom[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${star_qc_metrics[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${eset[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${f_data[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${p_data[@]}"; do
  cp -rL "$i" "$par_output/"
done

cp -rL "$par_html_report" "$par_output/"

cp -rL "$par_run_params" "$par_output/"
