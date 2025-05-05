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
# Keep track of all of the header files in this directory
#
set(_header_files
    ad_drivers.h
    ad_helpers.h
    ad_number_traits.h
    ad_number_types.h
    adolc_math.h
    adolc_number_types.h
    adolc_product_types.h
    sacado_math.h
    sacado_number_types.h
    sacado_product_types.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
