set(SRC_DOFS
  block_info.cc
  dof_accessor.cc
  dof_faces.cc
  dof_handler.cc
  dof_handler_policy.cc
  dof_levels.cc
  dof_objects.cc
  dof_renumbering.cc
  dof_tools.cc
  number_cache.cc
  )

add_library(OBJ_DOFS OBJECT ${SRC_DOFS})
