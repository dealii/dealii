#!/bin/bash
set -u

DEAL_II_SOURCE_DIR="${1}"
shift
CMAKE_CURRENT_SOURCE_DIR="${1}"
shift

#
# update .diff files from current .cc files
#

# Find all *.diff files in the current directory
diff_files="$(find ${CMAKE_CURRENT_SOURCE_DIR} -maxdepth 1 -type f -name "*.diff")"

for file in ${diff_files}; do
  # Extract the filename without the extension
  filename="$(basename "${file}" ".diff")"
  source_file="${DEAL_II_SOURCE_DIR}/examples/${filename}/${filename}.cc"
  modified_file="${filename}.cc"

  # Check if the corresponding .cc file exists
  if [[ -f "${source_file}" && -f "${modified_file}" ]]; then
    echo diff "${source_file}" "${modified_file}" "${file}"
    diff -c "${source_file}" "${modified_file}" > "${file}"
  else
    echo "No matching .cc files found for ${file}"
    exit 1
  fi
done

exit 0
