## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  mg_base.cc
  mg_dof_handler.cc
  mg_tools.cc
  mg_transfer_block.cc
  mg_transfer_component.cc
  mg_transfer_prebuilt.cc
  multigrid.cc
  )

SET(_inst
  mg_base.inst.in
  mg_dof_handler.inst.in
  mg_tools.inst.in
  mg_transfer_block.inst.in
  mg_transfer_component.inst.in
  mg_transfer_prebuilt.inst.in
  multigrid.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/multigrid/*.h
  )

DEAL_II_ADD_LIBRARY(obj_multigrid OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_multigrid "${_inst}")
