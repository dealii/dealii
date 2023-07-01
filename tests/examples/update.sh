#!/bin/bash

# use .diff files to create .cc files

# Find all *.diff files in the current directory
diff_files=$(find . -maxdepth 1 -type f -name "*.diff")

for file in $diff_files; do
  # Extract the filename without the extension
  filename=$(basename "$file" .diff)
  src="../../examples/$filename/$filename.cc"
  
  # Check if the corresponding .cc file exists
  if [[ -f "$src" ]]; then
    # Apply the diff and save it to a new file
    patch "$src" < "$file" -o "$filename.cc"
  else
      echo "No matching .cc file found for $diff_file"
      exit -1
  fi
done
