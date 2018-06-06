## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
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
  )

SET(_inst
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/optimization/rol/*.h
  )

DEAL_II_ADD_LIBRARY(obj_rol OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_rol "${_inst}")
