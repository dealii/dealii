## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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

#
# A small wrapper around ADD_LIBRARY that will define a target for each
# build type specified in DEAL_II_BUILD_TYPES. Only compatible with object
# targets (as used in the build system).
#
# It is assumed that the desired compilation configuration is set via
#   DEAL_II_CXX_FLAGS_${build}
#   DEAL_II_DEFINITIONS_${build}
#
# as well as the global (for all build types)
#   DEAL_II_CXX_FLAGS
#   DEAL_II_DEFINITIONS
#

macro(deal_ii_add_library _library)

  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)

    add_library(${_library}_${_build_lowercase}
      ${ARGN}
      )

    set_target_properties(${_library}_${_build_lowercase} PROPERTIES
      LINKER_LANGUAGE "CXX"
      )

    if(CMAKE_VERSION VERSION_LESS 3.9 OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

      set_target_properties(${_library}_${_build_lowercase} PROPERTIES
        COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
        COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
        )

    else()

      set(_flags "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
      separate_arguments(_flags)
      target_compile_options(${_library}_${_build_lowercase} PUBLIC
        $<$<COMPILE_LANGUAGE:CXX>:${_flags}>
        )

      target_compile_definitions(${_library}_${_build_lowercase}
        PUBLIC ${DEAL_II_DEFINITIONS} ${DEAL_II_DEFINITIONS_${_build}}
        )

    endif()

    set_property(GLOBAL APPEND PROPERTY DEAL_II_OBJECTS_${_build}
      "$<TARGET_OBJECTS:${_library}_${_build_lowercase}>"
      )
  endforeach()

endmacro()
