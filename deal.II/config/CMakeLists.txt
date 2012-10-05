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


CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )


###########################################################################
#                                                                         #
#               Setup and install the convenience macros:                 #
#                                                                         #
###########################################################################

LIST(APPEND DEAL_II_MACROS
  ${CMAKE_INSTALL_PREFIX}/${DEAL_II_CMAKE_MACROS_RELDIR}/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_INSTALL_PREFIX}/${DEAL_II_CMAKE_MACROS_RELDIR}/macro_deal_ii_setup_target.cmake
  )

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  DESTINATION ${DEAL_II_CMAKE_MACROS_RELDIR}
  COMPONENT library
  )


###########################################################################
#                                                                         #
#         Configure and install the cmake project configuration:          #
#                                                                         #
###########################################################################

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )

#
# Append target information:
#
FOREACH(build ${DEAL_II_BUILD_TYPES})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    "SET(DEAL_II_TARGET_${build} \"${DEAL_II_BASE_NAME}${DEAL_II_${build}_SUFFIX}\")\n"
    )
ENDFOREACH()

#
# Append feature configuration:
#
FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  "\n\n#\n# Feature configuration:\n#\n\n"
  )
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
      "SET(${var} ${${var}})\n"
      )
  ENDIF()
ENDFOREACH()

#
# Include convenience macros directly at the end of our project config:
#
FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  "\n\n#\n# Include some convenience macros directly in the project config:\n#\n\n"
  )
FOREACH(var ${DEAL_II_MACROS})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    "INCLUDE(${var})\n"
    )
ENDFOREACH()


CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

