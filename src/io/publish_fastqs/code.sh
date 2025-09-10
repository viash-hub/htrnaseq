#!/bin/bash

echo "Publishing $par_input -> $par_output"

echo
echo "Creating directory if it does not exist:"
mkdir -p "$par_output" && echo "$par_output created"

echo
echo "Copying files..."
IFS=";" read -ra input <<<$par_input

for i in "${input[@]}"; do
  cp -a "$i" "$par_output/"
done