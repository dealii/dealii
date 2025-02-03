## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_src
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/integrators/*.h
  )

define_object_library(object_integrators OBJECT ${_src} ${_header})


# Keep track of a complete list of all source files, based on the
# variables we have set above.  Prepend each source name by the
# current directory.
set(_source_files)
list(APPEND _source_files ${_src})
list(TRANSFORM _source_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these source files to the global list and propagate
# things back up one directory level
list(APPEND DEAL_II_SOURCE_FILES ${_source_files})
set(DEAL_II_SOURCE_FILES ${DEAL_II_SOURCE_FILES} PARENT_SCOPE)
