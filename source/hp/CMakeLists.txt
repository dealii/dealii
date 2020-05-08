## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_unity_include_src
  dof_faces.cc
  dof_level.cc
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  refinement.cc
  )

SET(_separate_src
  dof_handler.cc
  )

# concatenate all unity inclusion files in one file
SET(_n_includes_per_unity_file 15)

SETUP_SOURCE_LIST("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

SET(_inst
  dof_handler.inst.in
  fe_collection.inst.in
  fe_values.inst.in
  mapping_collection.inst.in
  refinement.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/hp/*.h
  )

DEAL_II_ADD_LIBRARY(obj_hp OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_hp "${_inst}")
