## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


#
# Recurse into the sub-directories:
#
add_subdirectory(deal.II)


#
# Sanity check: Make sure that the CMakeLists.txt files in the subdirectories
# really did list all header files.
#
file(GLOB_RECURSE _header_files
     LIST_DIRECTORIES false
     "${CMAKE_CURRENT_SOURCE_DIR}/deal.II/*.h")
foreach(_h ${DEAL_II_HEADER_FILES})
  list(REMOVE_ITEM _header_files ${_h})
endforeach()
list(LENGTH _header_files _count)

if (NOT ${_count} EQUAL 0)
  message(FATAL_ERROR
          "It looks like the CMakeLists.txt files in the subdirectories of\n"
          "include/ forgot to list the following ${_count} header files:\n"
          "  ${_header_files}")
endif()

# We will eventually need the names of header files to build interface
# module partitions. The subdirectories all reported the names of
# header files as global paths, so let us now strip the global parts
# so that all header files are reported relative to the 'include/'
# directory. That is, we need to strip CMAKE_CURRENT_SOURCE_DIR from
# all of the file names.
list(TRANSFORM DEAL_II_HEADER_FILES
     REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "")

# Propagate the list of header files back up one directory level
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)



#
# Add a rule for how to install the header files that are part of the source tree...
#
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deal.II
  DESTINATION ${DEAL_II_INCLUDE_RELDIR}
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  )

#
# ...and don't forget to install all generated header files, too:
#
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_INCLUDE_RELDIR}
  COMPONENT library
  FILES_MATCHING PATTERN "*.h"
  )
