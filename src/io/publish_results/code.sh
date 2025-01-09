#!/bin/bash

echo "Publishing $par_input -> $par_output"

echo
echo "Creating directory if it does not exist:"
mkdir -p "$par_output" && echo "$par_output created"

echo
echo "Copying files..."
IFS=";" read -ra star_output <<<$par_star_output

for i in "${star_output[@]}"; do
  cp -rL "$i" "$par_output/"
done

cp -rL "$par_nrReadsNrGenesPerChrom" "$par_output/"

cp -rL "$par_star_qc_metrics" "$par_output/"

cp -rL "$par_eset" "$par_output/"

cp -rL "$par_f_data" "$par_output/"

cp -rL "$par_p_data" "$par_output/"

cp -rL "$par_html_report" "$par_output/"
