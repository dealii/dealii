#!/bin/bash
#
# A script that can be used to create/update the .diff files in the
# current directory based on the original tutorial code and a modified
# version stored in a separate directory (typically in
# $BUILD/tests/examples/source). Call this script as in
#   ./update_diffs.sh /path/to/deal.II /path/to/deal.II/tests/examples
# This script is automatically executed when you say
#   make update_diffs
# in the build directory after setting up the tests.
#
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
    echo diff "${source_file}" "${modified_file}" \> "${file}"
    diff "${source_file}" "${modified_file}" > "${file}"
  else
    echo "No matching .cc files found for ${file}: source=${source_file} modified=${modified_file}"
  fi
done

exit 0
