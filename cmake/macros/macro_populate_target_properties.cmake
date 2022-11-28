## ---------------------------------------------------------------------
##
## Copyright (C) 2022 - 2022 by the deal.II authors
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
# populate_target_properties(<target> <build>)
#
# This function populate target properties according to (globally) defined
# DEAL_II_* variables. Specifically:
#
#   DEAL_II_LIBRARIES DEAL_II_LIBRARIES_<build>
#   DEAL_II_TARGETS DEAL_II_TARGETS_<build>
#     - populating the LINK_LIBRARIES target property
#   DEAL_II_INCLUDE_DIRS
#     - populating the INCLUDE_DIRECTORIES target property
#   DEAL_II_DEFINITIONS DEAL_II_DEFINITIONS_<build>
#     - populating the COMPILE_DEFINITIONS target property
#   DEAL_II_CXX_FLAGS DEAL_II_CXX_FLAGS_<build>
#     - populating the COMPILE_OPTIONS target property
#   DEAL_II_LINKER_FLAGS DEAL_II_LINKER_FLAGS_<build>
#     - populating the LINK_OPTIONS target property
#
# In addition the macro sets up a BUILD_INTERFACE and INSTALL_INTERFACE
# includes
#

function(populate_target_properties _target _build)

  set_target_properties(${_target} PROPERTIES LINKER_LANGUAGE "CXX")

  target_link_libraries(${_target}
    PUBLIC ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${_build}}
           ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_${_build}}
    )

  #
  # Include the current source directory of any target as private include
  # and add the contents of ${DEAL_II_INCLUDE_DIRS} as a public interface.
  #
  target_include_directories(${_target} BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${_target} SYSTEM PUBLIC ${DEAL_II_INCLUDE_DIRS})

  # Build directory specific includes:

  target_include_directories(${_target} PRIVATE
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${DEAL_II_BUNDLED_INCLUDE_DIRS}
    )

  # Interface includes:

  set(_includes
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    "$<INSTALL_INTERFACE:\${DEAL_II_INCLUDE_RELDIR}>"
    )
  foreach(_include ${DEAL_II_BUNDLED_INCLUDE_DIRS})
    list(APPEND _includes $<BUILD_INTERFACE:${_include}>)
  endforeach()
  if(NOT "${DEAL_II_BUNDLED_INCLUDE_DIRS}" STREQUAL "")
    list(APPEND _includes "$<INSTALL_INTERFACE:\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled>")
  endif()

  target_include_directories(${_target} INTERFACE ${_includes})

  target_compile_definitions(${_target}
    PUBLIC ${DEAL_II_DEFINITIONS} ${DEAL_II_DEFINITIONS_${_build}}
    )

  separate_arguments(_compile_options UNIX_COMMAND
    "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )
  target_compile_options(${_target}
    PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${_compile_options}>
    )

  separate_arguments(_link_options UNIX_COMMAND
    "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    )
  target_link_options(${_target} INTERFACE ${_link_options})

endfunction()
