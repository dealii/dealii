## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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
mark_as_advanced(file_cmd)

#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
foreach(_flags ${DEAL_II_USED_FLAGS})
  # Strip leading and trailing whitespace:
  string(STRIP "${${_flags}} ${${_flags}_SAVED}" ${_flags})
endforeach()

#
# Sanity check: The variables defined in DEAL_II_REMOVED_FLAGS must not be
# used during the configuration stage:
#
foreach(_flag ${DEAL_II_REMOVED_FLAGS})
  if(NOT "${${_flag}}" STREQUAL "")
    message(FATAL_ERROR
      "\nInternal configuration error: The variable ${_flag} was set to a "
      "non empty value during the configuration! (The corresponding "
      "DEAL_II_* variable should have been used.)\n"
      "${_flag}=\"${${_flag}}\"\n"
      )
  endif()
endforeach()

#
# Sanity check: Can we compile with the final setup?
#

foreach(build ${DEAL_II_BUILD_TYPES})
  macro(_check_linker_flags)
    check_compiler_setup(
      "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}}"
      "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
      DEAL_II_HAVE_USABLE_FLAGS_${build}
      ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${build}}
      ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_${build}}
      )
  endmacro()

  macro(_set_cache_variable _variable _value)
    set(${_variable} ${_value} CACHE INTERNAL "" FORCE)
    set(${_variable} ${_value})
  endmacro()

  macro(_drop_linker_flag _linker_flag _replacement_flag _variable)
    message(STATUS
      "Unable to compile a simple test program. "
      "Trying to drop \"${_linker_flag}\" from the linker flags."
      )
    foreach(_flags
        DEAL_II_LINKER_FLAGS DEAL_II_LINKER_FLAGS_${build}
        BASE_LINKER_FLAGS BASE_LINKER_FLAGS_${build}
        )
      string(REPLACE "${_linker_flag}" "${_replacement_flag}"
        ${_flags} "${${_flags}}"
        )
    endforeach()
    _set_cache_variable(_variable FALSE)
  endmacro()

  _check_linker_flags()

  if(NOT DEAL_II_HAVE_USABLE_FLAGS_${build} AND DEAL_II_COMPILER_HAS_FUSE_LD_MOLD)
    set(_replacement "")
    if(DEAL_II_COMPILER_HAS_FUSE_LD_LLD)
      set(_replacement "-fuse-ld=lld")
    elseif(DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
      set(_replacement "-fuse-ld=gold")
    endif()
    _drop_linker_flag(
      "-fuse-ld=mold" ${_replacement}
      DEAL_II_COMPILER_HAS_FUSE_LD_MOLD
      )
    _check_linker_flags()
  endif()

  if(NOT DEAL_II_HAVE_USABLE_FLAGS_${build} AND DEAL_II_COMPILER_HAS_FUSE_LD_LLD)
    set(_replacement "")
    if(DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
      set(_replacement "-fuse-ld=gold")
    endif()
    _drop_linker_flag(
      "-fuse-ld=lld" ${_replacement}
      DEAL_II_COMPILER_HAS_FUSE_LD_LLD
      )
    _check_linker_flags()
  endif()

  if(NOT DEAL_II_HAVE_USABLE_FLAGS_${build} AND DEAL_II_COMPILER_HAS_FUSE_LD_GOLD)
    _drop_linker_flag(
      "-fuse-ld=gold" ""
      DEAL_II_COMPILER_HAS_FUSE_LD_GOLD
      )
    _check_linker_flags()
  endif()

  if(NOT DEAL_II_HAVE_USABLE_FLAGS_${build})
    message(FATAL_ERROR "
  Configuration error: Cannot compile a test program with the final set of
  compiler and linker flags:
    CXX flags (${build}): ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}}
    LD flags  (${build}): ${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}
    LIBRARIES (${build}): ${DEAL_II_LIBRARIES};${DEAL_II_LIBRARIES_${build}}
  \n\n"
      )
  endif()
endforeach()

#
# Clean up deal.IITargets.cmake in the build directory:
#

file(REMOVE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  )


########################################################################
#                                                                      #
#                Query for git repository information:                 #
#                                                                      #
########################################################################

deal_ii_query_git_information("DEAL_II")

file(WRITE ${CMAKE_BINARY_DIR}/revision.log
"###
#
#  Git information:
#        Branch:    ${DEAL_II_GIT_BRANCH}
#        Revision:  ${DEAL_II_GIT_REVISION}
#        Timestamp: ${DEAL_II_GIT_TIMESTAMP}
#
###"
  )

