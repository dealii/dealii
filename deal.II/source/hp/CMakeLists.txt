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

ADD_LIBRARY(obj_hp OBJECT ${src})
macro_expand_instantiations(obj_hp "${inst_in_files}")
ADD_DEPENDENCIES(obj_hp deal_ii_target_dependencies)
