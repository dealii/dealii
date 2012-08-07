set(SRC_HP
  dof_faces.cc
  dof_handler.cc
  dof_levels.cc
  dof_objects.cc
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  )

add_library(OBJ_HP OBJECT ${SRC_HP})
