
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_DEFINITIONS name)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)

    GET_TARGET_PROPERTY(macro_definitions ${name}.${build_lowercase} COMPILE_DEFINITIONS)
    SET_TARGET_PROPERTIES(${name}.${build_lowercase} PROPERTIES
      COMPILE_DEFINITIONS "${ARGN};${macro_definitions}"
      )
  ENDFOREACH()

ENDMACRO()

