set(SRC_MULTIGRID
  mg_base.cc
  mg_dof_accessor.cc
  mg_dof_handler.cc
  mg_tools.cc
  mg_transfer_block.cc
  mg_transfer_component.cc
  mg_transfer_prebuilt.cc
  multigrid.cc
  )

add_library(OBJ_MULTIGRID OBJECT ${SRC_MULTIGRID})
