## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2025 by the deal.II authors
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
  mesh_worker.cc
  mesh_worker_info.cc
  mesh_worker_vector_selector.cc
  scratch_data.cc
  )

set(_inst
  mesh_worker_info.inst.in
  mesh_worker_vector_selector.inst.in
  scratch_data.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/meshworker/*.h
  )

define_object_library(object_meshworker OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_meshworker "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
