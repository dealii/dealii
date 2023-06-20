## ---------------------------------------------------------------------
##
## Copyright (C) 2022 - 2023 by the deal.II authors
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

  if(NOT "${_target}" MATCHES "^(object|bundled|${DEAL_II_TARGET_NAME})_")
    message(FATAL_ERROR
      "Internal error: The specified target name must begin with object_, "
      "bundled_, or ${DEAL_II_TARGET_NAME}_. Encountered: ${_target}"
      )
  endif()

  set(_visibility PRIVATE)
  if("${_target}" MATCHES "^${DEAL_II_TARGET_NAME}")
    set(_visibility PUBLIC)
  endif()

  set_target_properties(${_target} PROPERTIES LINKER_LANGUAGE "CXX")

  #
  # Add the contents of ${DEAL_II_INCLUDE_DIRS} as a public interface.
  #
  target_include_directories(${_target} SYSTEM ${_visibility} ${DEAL_II_INCLUDE_DIRS})

  # Build-directory specific includes:

  target_include_directories(${_target} PRIVATE
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    )
  target_include_directories(${_target} SYSTEM PRIVATE
    ${DEAL_II_BUNDLED_INCLUDE_DIRS}
    )

  # Interface includes:

  if("${_visibility}" STREQUAL "PUBLIC")
    set(_includes
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
      "$<INSTALL_INTERFACE:${DEAL_II_INCLUDE_RELDIR}>"
      )
    foreach(_include ${DEAL_II_BUNDLED_INCLUDE_DIRS})
      list(APPEND _includes $<BUILD_INTERFACE:${_include}>)
    endforeach()
    if(NOT "${DEAL_II_BUNDLED_INCLUDE_DIRS}" STREQUAL "")
      list(APPEND _includes "$<INSTALL_INTERFACE:${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled>")
    endif()
    target_include_directories(${_target} SYSTEM INTERFACE ${_includes})
  endif()

  target_compile_definitions(${_target} ${_visibility}
    ${DEAL_II_DEFINITIONS} ${DEAL_II_DEFINITIONS_${_build}}
    )

  #
  # Add target properties:
  #
  #  - set POSITION_INDEPENDENT_CODE to true to compile everything with
  #    the -fpic/-fPIC compiler flag. This ensures that we can link all
  #    object targets into a relocatable library at the end.
  #

  set_target_properties(${_target} PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    )

  #
  # Add compile and link options with private scope, and add the link
  # interface:
  #

  target_compile_flags(${_target} PRIVATE
    "${DEAL_II_WARNING_FLAGS} ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )

  get_property(_type TARGET ${_target} PROPERTY TYPE)
  if(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    target_link_flags(${_target} PRIVATE
      "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
      )
  endif()

  target_link_libraries(${_target} ${_visibility}
    ${DEAL_II_LIBRARIES} ${DEAL_II_LIBRARIES_${_build}}
    ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_${_build}}
    )
endfunction()
