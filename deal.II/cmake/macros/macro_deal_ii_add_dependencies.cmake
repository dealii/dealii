
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_DEPENDENCIES name target)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)
    ADD_DEPENDENCIES(${name}.${build_lowercase}
      ${target}.${build_lowercase}
      )
  ENDFOREACH()

ENDMACRO()
