## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2019 by the deal.II authors
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

include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

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

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/meshworker/*.h
  )

deal_ii_add_library(obj_meshworker OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_meshworker "${_inst}")
