## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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
# used during the configuration stage:
#
FOREACH(_flag ${DEAL_II_REMOVED_FLAGS})
  IF(NOT "${${_flag}}" STREQUAL "")
    MESSAGE(FATAL_ERROR
      "\nInternal configuration error: The variable ${_flag} was set to a "
      "non empty value during the configuration! (The corresponding "
      "DEAL_II_* variable should have been used.)\n"
      "${_flag}=\"${${_flag}}\"\n"
      )
  ENDIF()
ENDFOREACH()

#
# Save base configuration into variables BASE_* for later use in
# setup_write_config.cmake:
#
FOREACH(_suffix ${DEAL_II_STRING_SUFFIXES} ${DEAL_II_LIST_SUFFIXES})
  SET(BASE_${_suffix} ${DEAL_II_${_suffix}})
ENDFOREACH()

#
# Register features:
#
FOREACH(_feature ${DEAL_II_FEATURES})
  IF(DEAL_II_WITH_${_feature})
    FILTER_SYSTEM_LIBRARIES(${_feature})
    REGISTER_FEATURE(${_feature})
  ENDIF()
ENDFOREACH()

#
# Deduplicate entries one more time :-]
#
FOREACH(_suffix ${DEAL_II_LIST_SUFFIXES})
  IF(_suffix MATCHES "INCLUDE_DIRS$")
    REMOVE_DUPLICATES(DEAL_II_${_suffix})
  ELSE()
    REMOVE_DUPLICATES(DEAL_II_${_suffix} REVERSE)
  ENDIF()
ENDFOREACH()

#
# Sanity check: Can we compile with the final setup?
#

FOREACH(build ${DEAL_II_BUILD_TYPES})
# FIXME: until https://github.com/dealii/dealii/issues/3686 is resolved,
# the simple tests below hangs and renders deal.II unusable.
  SET(DEAL_II_HAVE_USABLE_FLAGS_${build} TRUE)
#  CHECK_COMPILER_SETUP(
#    "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}}"
#    "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
#    DEAL_II_HAVE_USABLE_FLAGS_${build}
#    ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${build}}
#    )

  IF(NOT DEAL_II_HAVE_USABLE_FLAGS_${build})
    MESSAGE(FATAL_ERROR "
  Configuration error: Cannot compile a test program with the final set of
  compiler and linker flags:
    CXX flags (${build}): ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}}
    LD flags  (${build}): ${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}
    LIBRARIES (${build}): ${DEAL_II_LIBRARIES};${DEAL_II_LIBRARIES_${build}}
  \n\n"
      )
  ENDIF()
ENDFOREACH()

#
# Clean up deal.IITargets.cmake in the build directory:
#
FILE(REMOVE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  )
