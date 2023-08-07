#!/bin/bash

directory="./data/CSVs/"

# Check if the directory exists
if [ ! -d "$directory" ]; then
  echo "Error: Directory not found."
  exit 1
fi

extract_occs() {
  local file="$1"
  echo "Processing file: $file"
  Rscript preprocessing.R "$file"
}

export -f extract_occs

find "$directory" -type f -name "*.csv" | xargs -I {} -P 6 bash -c 'extract_occs "$@"' _ {}
