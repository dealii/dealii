## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
#                      Finalize the configuration:                     #
#                                                                      #
########################################################################

#
# Hide some cmake specific cached variables. This is annoying...
#
MARK_AS_ADVANCED(file_cmd)

#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
FOREACH(_flags ${DEAL_II_USED_FLAGS})
  # Strip leading and trailing whitespace:
  STRING(STRIP "${${_flags}} ${${_flags}_SAVED}" ${_flags})
ENDFOREACH()

#
# Sanity check: The variables defined in DEAL_II_REMOVED_FLAGS must not be
# used during the comfiguration stage:
#
FOREACH(_flag ${DEAL_II_REMOVED_FLAGS})
  IF(NOT "${_flag}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "\nInternal configuration error: The variable ${_flag} was set to a "
      "non empty value during the configuration! (The corresponding "
      "DEAL_II_* variable should have been used.)\n"
      "${_flag}=\"${${_flag}}\"\n"
      )
  ENDIF()
ENDFOREACH()

#
# Deduplicate entries in DEAL_II_USER_INCLUDE_DIRS
#
REMOVE_DUPLICATES(DEAL_II_USER_INCLUDE_DIRS)

#
# Deduplicate entries in DEAL_II_EXTERNAL_LIBRARIES(_...)
# in reverse order:
#
REMOVE_DUPLICATES(DEAL_II_EXTERNAL_LIBRARIES REVERSE)

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  REMOVE_DUPLICATES(DEAL_II_EXTERNAL_LIBRARIES_${_build} REVERSE)
ENDFOREACH()

#
# Cleanup some files used for storing the names of all object targets that
# will be bundled to the deal.II library.
# (Right now, i.e. cmake 2.8.8, this is the only reliable way to get
# information into a global scope...)
#
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER "${_build}" _build_lowercase)
  FILE(REMOVE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${_build_lowercase}
    )
ENDFOREACH()
FILE(WRITE
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_source_includes
  ""
  )

#
# Cleanup deal.IITargets.cmake in the build directory:
#
FILE(WRITE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  ""
  )
