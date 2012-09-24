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
DEAL_II_ADD_DEPENDENCIES(obj_multigrid deal_ii_target_dependencies)
