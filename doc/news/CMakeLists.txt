## -----------------------------------------------------------------------------
##
## Copyright (C) 2016 - 2022 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------


#
# Set up all necessary bits for the changelog generation
#
if(DEAL_II_COMPONENT_DOCUMENTATION)

  message(STATUS "")
  message(STATUS "Setting up changelog generation")

  #
  # Generate changes.h.
  #
  add_custom_command(
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

  add_custom_target(changelog ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/changes.h
  )

  add_dependencies(documentation changelog)

endif()
