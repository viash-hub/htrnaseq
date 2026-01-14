#!/bin/bash

echo "Publishing $par_input -> $par_output"

echo
echo "Creating directory if it does not exist:"
mkdir -p "$par_output" && echo "$par_output created"

echo
if [ "$par_skip" == "true" ]; then
  echo "Skipping publishing of fastq files."
  # Do not copy
  # But we must ensure the output parameter is satisfied if it's required by Viash/Nextflow
  # Since par_output default is "$id/", and direction is output, Viash expects this directory to be created.
  # We will create an empty directory.
  mkdir -p "$par_output"
else 
  echo "Copying files..."
  IFS=";" read -ra input <<<$par_input

  for i in "${input[@]}"; do
    cp -a --keep-directory-symlink "$i" "$par_output/"
  done
fi