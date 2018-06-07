## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

# Auxiliary functions
FUNCTION(CAT IN_FILE OUT_FILE INDENT)
  FILE(READ ${IN_FILE} CONTENTS)
  IF(${INDENT} MATCHES "TRUE")
    FILE(STRINGS ${IN_FILE} LINESTMP)
    FOREACH(LINETMP ${LINESTMP})
      FILE(APPEND ${OUT_FILE} "  ${LINETMP}\n")
    ENDFOREACH()
  ELSE()
    FILE(APPEND ${OUT_FILE} "${CONTENTS}")
  ENDIF()
ENDFUNCTION()

FUNCTION(PROCESS IN_DIR OUT_FILE)
  FILE(APPEND ${OUT_FILE} "<ol>\n")
  FILE(GLOB ENTRY_LIST ${IN_DIR}/[0-9]*)
  LIST(SORT ENTRY_LIST)
  LIST(REVERSE ENTRY_LIST)
  FOREACH(ENTRY ${ENTRY_LIST})
    FILE(APPEND ${OUT_FILE} "\n <li>\n")
    CAT(${ENTRY} ${OUT_FILE} "TRUE")
    FILE(APPEND ${OUT_FILE} " </li>\n")
  ENDFOREACH()
  FILE(APPEND ${OUT_FILE} "\n</ol>\n")
ENDFUNCTION()

# Generate 'changes.h'.

# First, create a file 'changes.h.in' based on all changelog fragments.
SET(OUTPUT_FILE_TEMP "${OUTPUT_FILE}.in")
FILE(WRITE ${OUTPUT_FILE_TEMP} "")
CAT    (${CMAKE_CURRENT_SOURCE_DIR}/header
        ${OUTPUT_FILE_TEMP} "FALSE")
CAT    (${CMAKE_CURRENT_SOURCE_DIR}/header_incompatibilities
        ${OUTPUT_FILE_TEMP} "FALSE")
PROCESS(${CMAKE_CURRENT_SOURCE_DIR}/incompatibilities
        ${OUTPUT_FILE_TEMP})
CAT    (${CMAKE_CURRENT_SOURCE_DIR}/header_major
        ${OUTPUT_FILE_TEMP} "FALSE")
PROCESS(${CMAKE_CURRENT_SOURCE_DIR}/major
        ${OUTPUT_FILE_TEMP})
CAT    (${CMAKE_CURRENT_SOURCE_DIR}/header_minor
        ${OUTPUT_FILE_TEMP} "FALSE")
PROCESS(${CMAKE_CURRENT_SOURCE_DIR}/minor
        ${OUTPUT_FILE_TEMP})
CAT    (${CMAKE_CURRENT_SOURCE_DIR}/footer
        ${OUTPUT_FILE_TEMP} "FALSE")

# Copy it over to 'changes.h' but only touch the time stamp
# if the file actually changed (this is what CONFIGURE_FILE does).
MESSAGE(STATUS "Generating changelog file: ${OUTPUT_FILE}")
CONFIGURE_FILE(${OUTPUT_FILE_TEMP} ${OUTPUT_FILE} COPYONLY)
