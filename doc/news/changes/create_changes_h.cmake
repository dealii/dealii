## -----------------------------------------------------------------------------
##
## Copyright (C) 2016 - 2022 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

# Auxiliary functions
function(cat IN_FILE OUT_FILE INDENT)
  file(READ ${IN_FILE} CONTENTS)
  if(${INDENT} MATCHES "TRUE")
    file(STRINGS ${IN_FILE} LINESTMP ENCODING UTF-8)
    foreach(LINETMP ${LINESTMP})
      file(APPEND ${OUT_FILE} "  ${LINETMP}\n")
    endforeach()
  else()
    file(APPEND ${OUT_FILE} "${CONTENTS}")
  endif()
endfunction()

# Given a directory, take all of the changelog files in that
# directory (assumed to be starting with a date, i.e., number)
# and process them into an HTML-formatted list. Write that into
# an output file.
function(process IN_DIR OUT_FILE)
  file(APPEND ${OUT_FILE} "<ol>\n")

  # Get all of the files in this directory and sort the
  # list of files newest-to-oldest
  file(GLOB ENTRY_LIST ${IN_DIR}/[0-9]*)
  list(SORT ENTRY_LIST)
  list(REVERSE ENTRY_LIST)

  # Then concatenate as a bullet point list. Exclude entries
  # called 'dummy' that we place in each directory at the
  # beginning of each release process to ensure no directory
  # is empty. Also exclude files that are back-up files as
  # indicated by having a trailing tilde.
  list(LENGTH ENTRY_LIST n_entries)
  foreach(ENTRY ${ENTRY_LIST})
    if (((n_entries GREATER 1) AND (ENTRY MATCHES ".*dummy$"))
        OR (ENTRY MATCHES ".*~$"))
      continue()
    endif()

    file(APPEND ${OUT_FILE} "\n <li>\n")
    cat(${ENTRY} ${OUT_FILE} "TRUE")
    file(APPEND ${OUT_FILE} " </li>\n")
  endforeach()

  file(APPEND ${OUT_FILE} "\n</ol>\n")
endfunction()

# Generate 'changes.h'.

# First, create a file 'changes.h.in' based on all changelog fragments.
set(OUTPUT_FILE_TEMP "${OUTPUT_FILE}.in")
file(WRITE ${OUTPUT_FILE_TEMP} "")

# Then concatenate its various building blocks, where necessary
# creating these building blocks via the 'process()' function
# above.
cat    (${CMAKE_CURRENT_SOURCE_DIR}/header
        ${OUTPUT_FILE_TEMP} "FALSE")
cat    (${CMAKE_CURRENT_SOURCE_DIR}/header_incompatibilities
        ${OUTPUT_FILE_TEMP} "FALSE")
process(${CMAKE_CURRENT_SOURCE_DIR}/incompatibilities
        ${OUTPUT_FILE_TEMP})
cat    (${CMAKE_CURRENT_SOURCE_DIR}/header_major
        ${OUTPUT_FILE_TEMP} "FALSE")
process(${CMAKE_CURRENT_SOURCE_DIR}/major
        ${OUTPUT_FILE_TEMP})
cat    (${CMAKE_CURRENT_SOURCE_DIR}/header_minor
        ${OUTPUT_FILE_TEMP} "FALSE")
process(${CMAKE_CURRENT_SOURCE_DIR}/minor
        ${OUTPUT_FILE_TEMP})
cat    (${CMAKE_CURRENT_SOURCE_DIR}/footer
        ${OUTPUT_FILE_TEMP} "FALSE")

# Copy it over to 'changes.h' but only touch the time stamp
# if the file actually changed (this is what CONFIGURE_FILE does).
message(STATUS "Generating changelog file: ${OUTPUT_FILE}")
configure_file(${OUTPUT_FILE_TEMP} ${OUTPUT_FILE} COPYONLY)
