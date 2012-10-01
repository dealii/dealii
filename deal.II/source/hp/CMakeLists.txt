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
  dof_faces.cc
  dof_handler.cc
  dof_levels.cc
  dof_objects.cc
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  )

SET(inst_in_files
  dof_handler.inst.in
  fe_collection.inst.in
  fe_values.inst.in
  mapping_collection.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_hp OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_hp "${inst_in_files}")

