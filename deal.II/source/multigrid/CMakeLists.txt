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

SET(src
  mg_base.cc
  mg_dof_accessor.cc
  mg_dof_handler.cc
  mg_tools.cc
  mg_transfer_block.cc
  mg_transfer_component.cc
  mg_transfer_prebuilt.cc
  multigrid.cc
  )

SET(inst_in_files
  mg_dof_accessor.inst.in
  mg_dof_handler.inst.in
  mg_tools.inst.in
  mg_transfer_block.inst.in
  mg_transfer_component.inst.in
  mg_transfer_prebuilt.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_multigrid OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_multigrid "${inst_in_files}")

