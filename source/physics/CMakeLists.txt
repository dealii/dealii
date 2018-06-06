## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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

# Expand instantiations in subdirectories
ADD_SUBDIRECTORY("elasticity")

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  transformations.cc
  )

SET(_inst
  transformations.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/physics/*.h
  )

DEAL_II_ADD_LIBRARY(obj_physics OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_physics "${_inst}")
