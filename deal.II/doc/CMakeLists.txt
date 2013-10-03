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




#
# Set up all necessary bits for the documentation
#
IF(DEAL_II_COMPONENT_DOCUMENTATION)

  #
  # Required packages for the document generation:
  #
  FIND_PACKAGE(Perl)
  FIND_PACKAGE(Doxygen)

  #
  # Do we have all necessary dependencies?
  #
  IF(NOT PERL_FOUND)
    MESSAGE(FATAL_ERROR
      "Could not find a perl installation which is required for building the documentation"
      )
  ENDIF()
  IF(NOT DOXYGEN_FOUND OR NOT DOXYGEN_DOT_FOUND)
    MESSAGE(FATAL_ERROR
      "Could not find doxygen and dot which is required for building the documentation"
      )
  ENDIF()

  #
  # Generate the source documentation via doxygen:
  #
  ADD_SUBDIRECTORY(doxygen)

  #
  # A custom target with the same name as the component:
  #
  ADD_DEPENDENCIES(documentation doxygen)

  #
  # Install the static elements of the html documentation:
  #
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_DOCHTML_RELDIR}
    COMPONENT documentation
    PATTERN "CMakeLists.txt" EXCLUDE
    PATTERN "doxygen" EXCLUDE
    PATTERN "*.in" EXCLUDE
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/title.html.in
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    )
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    DESTINATION ${DEAL_II_DOCHTML_RELDIR}
    COMPONENT documentation
    )

ENDIF(DEAL_II_COMPONENT_DOCUMENTATION)

#
# Always install a minimalistic README and LICENSE file:
#

# Bloody hack to get an empty relative path:
IF(NOT "${DEAL_II_DOCREADME_RELDIR}" STREQUAL "")
  SET(_foo ${DEAL_II_DOCREADME_RELDIR})
ELSE()
  SET(_foo \${__empty})
ENDIF()

SET(_destination)
INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/README
  ${CMAKE_SOURCE_DIR}/LICENSE
  DESTINATION ${_foo}
  COMPONENT library
  )
