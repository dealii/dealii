## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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

set(_unity_include_src
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  refinement.cc
  )

set(_separate_src
  )

# concatenate all unity inclusion files in one file
set(_n_includes_per_unity_file 15)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  fe_collection.inst.in
  fe_values.inst.in
  mapping_collection.inst.in
  refinement.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/hp/*.h
  )

deal_ii_add_library(obj_hp OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_hp "${_inst}")
