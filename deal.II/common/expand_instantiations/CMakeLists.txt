#TODO template-arguments.in

ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)

GET_TARGET_PROPERTY (
  expand_instantiations_exe
  expand_instantiations
  LOCATION
  )

MACRO(macro_expand_instantiations inst_in_files)
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
ENDMACRO()
