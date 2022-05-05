## ---------------------------------------------------------------------
##
## Copyright (C) 2022 by the deal.II authors
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


SET(_src
surface_mesh.cc
)



SET(_inst
surface_mesh.inst.in
)

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/cgal/*.h
  )

DEAL_II_ADD_LIBRARY(obj_cgal OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_cgal "${_inst}")
