## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_src)
set(_inst)

if (DEAL_II_WITH_ADOLC OR DEAL_II_TRILINOS_WITH_SACADO)
  set(_src
    ${_src}
    ad_drivers.cc
    ad_helpers.cc
    )

  set(_inst
    ${_inst}
    ad_drivers.inst.in
    ad_drivers.inst1.in
    ad_drivers.inst2.in
    ad_helpers.inst.in
    ad_helpers.inst1.in
    ad_helpers.inst2.in
    )
endif()

if(DEAL_II_WITH_ADOLC)
  set(_src
    ${_src}
    adolc_number_types.cc
    )

  set(_inst
    ${_inst}
    adolc_number_types.inst.in
    )
endif()


if(DEAL_II_TRILINOS_WITH_SACADO)
  set(_src
    ${_src}
    sacado_number_types.cc
    )

  set(_inst
    ${_inst}
    sacado_number_types.inst1.in
    sacado_number_types.inst2.in
    )
endif()


file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/differentiation/ad/*.h
  )

define_object_library(object_differentiation_ad OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_differentiation_ad "${_inst}")


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
