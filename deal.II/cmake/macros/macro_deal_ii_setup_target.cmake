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
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_SETUP_TARGET(target)
#
# This macro sets the necessary include directories, linker flags, compile
# definitions and the external libraries the target will be linked against.
#

MACRO(DEAL_II_SETUP_TARGET target)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDE)
    MESSAGE(FATAL_ERROR
      "DEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use."
      )
  ENDIF()

  CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

  SET_TARGET_PROPERTIES(${target} PROPERTIES
    INCLUDE_DIRECTORIES
      "${DEAL_II_EXTERNAL_INCLUDE_DIRS};${DEAL_II_INCLUDE_DIRS}"
    LINK_FLAGS
      "${DEAL_II_LINKER_FLAGS}"
    LINK_FLAGS_DEBUG
      "${DEAL_II_LINKER_FLAGS_DEBUG}"
    LINK_FLAGS_RELEASE
      "${DEAL_II_LINKER_FLAGS_RELEASE}"
    COMPILE_DEFINITIONS
      "${DEAL_II_USER_DEFINITIONS}"
    COMPILE_DEFINITIONS_DEBUG
      "${DEAL_II_USER_DEFINITIONS_DEBUG}"
    COMPILE_DEFINITIONS_RELEASE
      "${DEAL_II_USER_DEFINITIONS_RELEASE}"
    )

    TARGET_LINK_LIBRARIES(${target}
        ${DEAL_II_LIBRARIES}
        ${DEAL_II_EXTERNAL_LIBRARIES}
      )

ENDMACRO()

