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

SET(src
  block_info.cc
  dof_accessor.cc
  dof_faces.cc
  dof_handler.cc
  dof_handler_policy.cc
  dof_objects.cc
  dof_renumbering.cc
  dof_tools.cc
  number_cache.cc
  )

SET(inst_in_files
  block_info.inst.in
  dof_accessor.inst.in
  dof_handler.inst.in
  dof_handler_policy.inst.in
  dof_objects.inst.in
  dof_renumbering.inst.in
  dof_tools.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_dofs OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_dofs "${inst_in_files}")

