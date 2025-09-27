## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# collect_source_and_header_files("<source files>" "<header files>")
#
# This macro appends a list of sources and header files to the global
# properties DEAL_II_SOURCE_FILES and DEAL_II_HEADER_FILES.
#

function(collect_source_and_header_files _sources _headers)

  #
  # Drop empty strings and check whether the first list element is an
  # absolute path. If not, then prepend the CMAKE_CURRENT_SOURCE_DIR to the
  # path:
  #
  foreach(_list_name _sources _headers)
    list(REMOVE_ITEM ${_list_name} "")

    if(NOT "${${_list_name}}" STREQUAL "")
      list(GET ${_list_name} 0 _first_element)
      if(NOT IS_ABSOLUTE ${_first_element})
        list(TRANSFORM ${_list_name} PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
      endif()
    endif()
  endforeach()

  set_property(GLOBAL APPEND PROPERTY DEAL_II_SOURCE_FILES ${_sources})
  set_property(GLOBAL APPEND PROPERTY DEAL_II_HEADER_FILES ${_headers})
endfunction()
