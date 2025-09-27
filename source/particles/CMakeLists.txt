## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
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
  data_out.cc
  particle.cc
  particle_handler.cc
  generators.cc
  property_pool.cc
  utilities.cc
  )

set(_inst
  data_out.inst.in
  particle.inst.in
  particle_handler.inst.in
  generators.inst.in
  utilities.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/particles/*.h
  )

define_object_library(object_particles OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_particles "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
