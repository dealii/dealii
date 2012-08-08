SET(src_distributed
  grid_refinement.cc
  solution_transfer.cc
  tria.cc
  )

SET(inst_in_files
  grid_refinement.inst.in
  solution_transfer.inst.in
  tria.inst.in
  )

GET_TARGET_PROPERTY (
  expand_instantiations_exe
  expand_instantiations
  LOCATION
  )

FOREACH (inst_in_file ${inst_in_files})
  STRING(REGEX REPLACE "\\.in$" "" inst_file "${inst_in_file}" )
  ADD_CUSTOM_COMMAND (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
    DEPENDS
      expand_instantiations
    COMMAND ${expand_instantiations_exe}
    ARGS
      ~/template-arguments
      < ${CMAKE_CURRENT_SOURCE_DIR}/${inst_in_file}
      > ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
    )
ENDFOREACH()

ADD_LIBRARY(obj_distributed OBJECT ${src_distributed})

