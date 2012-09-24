#
# A macro for the inst.in file expansion
#
# Usage:
#     EXPAND_INSTANTATIONS(target inst_in_files)
#
# Options:
#
# target
#
#    a target that will depend on the generation of all .inst files.
#    (To ensure that all .inst files are generated prior to compiling.)
#
# inst_in_files
#
#    a list of inst.in files that will be expanded
#

MACRO(EXPAND_INSTANTIATIONS target inst_in_files)
  FOREACH (inst_in_file ${inst_in_files})
    STRING(REGEX REPLACE "\\.in$" "" inst_file "${inst_in_file}" )

    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
      DEPENDS expand_instantiations
              ${CMAKE_CURRENT_SOURCE_DIR}/${inst_in_file}
      COMMAND expand_instantiations
      ARGS ${CMAKE_BINARY_DIR}/config/template-arguments
           < ${CMAKE_CURRENT_SOURCE_DIR}/${inst_in_file}
           > ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
      )

    LIST(APPEND inst_targets ${CMAKE_CURRENT_BINARY_DIR}/${inst_file})
  ENDFOREACH()

  #
  # Define a custom target that depends on the generation of all inst.in
  # files.
  #
  ADD_CUSTOM_TARGET(${target}.inst ALL DEPENDS ${inst_targets})

  #
  # Add a dependency to target so that target.inst is fully generated
  # before target will be processed.
  #
  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    ADD_DEPENDENCIES(${target}.g ${target}.inst)
  ENDIF()
  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    ADD_DEPENDENCIES(${target} ${target}.inst)
  ENDIF()

ENDMACRO()

