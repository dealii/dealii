## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_unity_include_src
  utilities.cc
  )

# include all files in the unity file
set(_n_includes_per_unity_file 20)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  utilities.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/vtk/*.h
  )

define_object_library(object_vtk OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_vtk "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
