MACRO(macro_expand_instantiations target inst_in_files)

  FOREACH (inst_in_file ${inst_in_files})
    STRING(REGEX REPLACE "\\.in$" "" inst_file "${inst_in_file}" )

    ADD_CUSTOM_COMMAND (
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
      DEPENDS expand_instantiations
              ${CMAKE_CURRENT_SOURCE_DIR}/${inst_in_file}
      COMMAND expand_instantiations
      ARGS ${CMAKE_BINARY_DIR}/common/expand_instantiations/template-arguments
           < ${CMAKE_CURRENT_SOURCE_DIR}/${inst_in_file}
           > ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
      )

    SET(inst_targets
      ${CMAKE_CURRENT_BINARY_DIR}/${inst_file}
      ${inst_targets}
      )
  ENDFOREACH()

  ADD_CUSTOM_TARGET(${target}.inst ALL DEPENDS ${inst_targets})

  ADD_DEPENDENCIES(${target} ${target}.inst)

ENDMACRO()
