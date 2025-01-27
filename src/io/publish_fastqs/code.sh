#!/bin/bash

echo "Publishing $par_input -> $par_output"

echo
echo "Creating directory if it does not exist:"
mkdir -p "$par_output" && echo "$par_output created"

echo
echo "Copying files..."
IFS=";" read -ra input_r1 <<<$par_input_r1
IFS=";" read -ra input_r2 <<<$par_input_r2

for i in "${input_r1[@]}"; do
  cp -rL "$i" "$par_output/"
done

for i in "${input_r2[@]}"; do
  cp -rL "$i" "$par_output/"
done
