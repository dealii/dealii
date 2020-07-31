## ---------------------------------------------------------------------
##
## Copyright (C) 2020 by the deal.II authors
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
  fe_lib.cc
  grid_generator.cc
  polynomials.cc
  quadrature_lib.cc
  )

SET(_separate_src
  )

# include all files in the unity file
SET(_n_includes_per_unity_file 20)

SETUP_SOURCE_LIST("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

SET(_inst
  fe_lib.inst.in
  grid_generator.inst.in
)

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/tet/*.h
  )

DEAL_II_ADD_LIBRARY(obj_simplex OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_simplex "${_inst}")
