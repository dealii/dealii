#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# See macro_deal_ii_add_library.cmake for an explanation of this macro
#

MACRO(DEAL_II_ADD_C_LIBRARY _library)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    ADD_LIBRARY(${_library}.${_build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${_library}.${_build_lowercase} PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_${_build}}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_${_build}}"
      )

    FILE(APPEND
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${_build_lowercase}
      "$<TARGET_OBJECTS:${_library}.${_build_lowercase}>\n"
      )
  ENDFOREACH()

ENDMACRO()
