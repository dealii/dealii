#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/hp/*.h
  )

SET(_src
  dof_faces.cc
  dof_handler.cc
  dof_levels.cc
  dof_objects.cc
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  )

SET(_inst
  dof_handler.inst.in
  fe_collection.inst.in
  fe_values.inst.in
  mapping_collection.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_hp OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_hp "${_inst}")

