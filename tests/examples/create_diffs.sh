#!/bin/bash

# create .diff files from current .cc files

# Find all *.diff files in the current directory
diff_files=$(find . -maxdepth 1 -type f -name "*.diff")

for file in $diff_files; do
  # Extract the filename without the extension
  filename=$(basename "$file" .diff)
  src="../../examples/$filename/$filename.cc"

  # Check if the corresponding .cc file exists
  if [[ -f "$src" ]]; then
    echo "creating $file from $src and $filename.cc"
    diff "$src" "$filename.cc" > "$file"
  else
    echo "No matching .cc file found for $file"
    exit 1
  fi
done
