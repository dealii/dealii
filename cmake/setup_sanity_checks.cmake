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
#                            Sanity checks:                            #
#                                                                      #
########################################################################

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
