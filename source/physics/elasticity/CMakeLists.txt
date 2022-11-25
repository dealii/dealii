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


include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

set(_src
  kinematics.cc
  standard_tensors.cc
  )

set(_inst
  kinematics.inst.in
  standard_tensors.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/physics/elasticity/*.h
  )

deal_ii_add_library(obj_physics_elasticity OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_physics_elasticity "${_inst}")
