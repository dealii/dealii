## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2020 by the deal.II Authors
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
# Set up all necessary bits for the changelog generation
#
IF(DEAL_II_COMPONENT_DOCUMENTATION)

  MESSAGE(STATUS "")
  MESSAGE(STATUS "Setting up changelog generation")

  #
  # Generate changes.h.
  #
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/changes.h
    COMMAND ${CMAKE_COMMAND} -DOUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/changes.h
                             -P ${CMAKE_CURRENT_SOURCE_DIR}/changes/create_changes_h.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/changes
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/header
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/header_incompatibilities
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/header_major
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/header_minor
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/footer
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/major
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/minor
      ${CMAKE_CURRENT_SOURCE_DIR}/changes/incompatibilities
  )

  ADD_CUSTOM_TARGET(changelog ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/changes.h
  )

  ADD_DEPENDENCIES(documentation changelog)

ENDIF(DEAL_II_COMPONENT_DOCUMENTATION)

