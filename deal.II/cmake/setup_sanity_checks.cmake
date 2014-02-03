## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------


########################################################################
#                                                                      #
#                            Sanity checks:                            #
#                                                                      #
########################################################################

#
# A quick test whether we're able to successfully link with the given
# compiler and linker flags and the given library link interface:
#

MESSAGE(STATUS "")
MESSAGE(STATUS "Sanity checks.")

FOREACH(_build ${DEAL_II_BUILD_TYPES})

  FOREACH(_var
    CXX_FLAGS CXX_FLAGS_${_build}
    LINKER_FLAGS LINKER_FLAGS_${_build}
    LIBRARIES LIBRARIES_${_build}
    )
    IF(NOT "${DEAL_II_${_var}}" STREQUAL "${CACHED_DEAL_II_${_var}_${_build}}")
      UNSET(DEAL_II_SANITY_CHECK_${_build} CACHE)
      SET(CACHED_DEAL_II_${_var}_${_build} "${DEAL_II_${_var}}" CACHE INTERNAL "" FORCE)
    ENDIF()
  ENDFOREACH()

  RESET_CMAKE_REQUIRED()
  ADD_FLAGS(CMAKE_REQUIRED_FLAGS
    "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES
    "${DEAL_II_LINKER_FLAGS} ${DEAL_II_CXX_LINKER_${_build}}"
    )
  LIST(APPEND CMAKE_REQUIRED_LIBRARIES
    ${DEAL_II_LIBRARIES}
    ${DEAL_II_LIBRARIES_${_build}}
    )
  CHECK_CXX_SOURCE_COMPILES("int main(){ return 0; }"
    DEAL_II_SANITY_CHECK_${_build}
    )
  RESET_CMAKE_REQUIRED()

  IF(NOT DEAL_II_SANITY_CHECK_${_build})
    UNSET(DEAL_II_SANITY_CHECK_${_build} CACHE)
    MESSAGE(FATAL_ERROR "
  Configuration error: Cannot compile and link with the current set of
  compiler flags, linker flags and libraries!

  Please check the test output given at the end of
  CMakeFiles/CMakeError.log and consult detailed.log for the current
  configuration.\n\n"
    )
  ENDIF()
ENDFOREACH()
