#!/bin/bash

set -eo pipefail

echo "Publishing results to multiple output directories"

# Create output directories for multiple files
echo "Creating output directories..."

declare -A path_pars_dirs=(
  ["par_star_output_dir"]="par_star_output"
  ["par_nrReadsNrGenesPerChrom_dir"]="par_nrReadsNrGenesPerChrom"
  ["par_star_qc_metrics_dir"]="par_star_qc_metrics"
  ["par_eset_dir"]="par_eset"
  ["par_f_data_dir"]="par_f_data"
  ["par_p_data_dir"]="par_p_data"
)

declare -A path_pars_files=(
  ["par_html_report_output"]="par_html_report"
  ["par_run_params_output"]="par_run_params"
)

echo "Canonicalizing output paths."
all_path_args=( "${!path_pars_files[@]}" "${!path_pars_dirs[@]}" )
for par in ${all_path_args[@]}; do
    curr_val="${!par}"
    printf "\t%s\n" "Canonicalizing path '$curr_val'"
    new_value=$(realpath --canonicalize-missing "$curr_val")
    printf "\t%s\n" "New output path for '$par': '$new_value'"
    declare -g "$par=$new_value"
done

echo "Creating output directories"
for par in ${!path_pars_dirs[@]}; do
    curr_val="${!par}"
    mkdir -p "$curr_val" && printf "\t%s\n" "$curr_val created."
done
for par in ${!path_pars_files[@]}; do
    curr_val="${!par}"
    file_dir=$(dirname "$curr_val")
    mkdir -p "$file_dir" && printf "\t%s\n" "$file_dir created."
done

echo "Copying output for inputs that are directories"
for par in ${!path_pars_dirs[@]}; do

  if [ "$par_skip" == "true" ] && [[ "$par" == "par_star_output_dir" || "$par" == "par_eset_dir" || "$par" == "par_f_data_dir" || "$par" == "par_p_data_dir" ]]; then
    printf "\t%s\n" "Skipping output for '$par' as --skip is enabled."
    continue
  fi

  printf "\t%s\n" "Copying output for '$par'"
  output_val="${!par}"
  input_par_name="${path_pars_dirs[${par}]}"
  input_val="${!input_par_name}"
  IFS=";" read -ra output_files_list <<<$input_val
  for i in "${output_files_list[@]}"; do 
      printf "\t%s\n" "Copying '$i' to '$output_val'" 
      cp -a --keep-directory-symlink "$i" "$output_val/"
  done
done

echo "Copying single files directly..."
for par in ${!path_pars_files[@]}; do
  output_val="${!par}"
  input_par_name="${path_pars_files[${par}]}"
  input_val="${!input_par_name}"
  cp -a --keep-directory-symlink "$input_val" "$output_val"
done

echo "Publishing completed successfully!"
