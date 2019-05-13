## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II Authors
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
# Set up all necessary bits for the documentation
#
IF(DEAL_II_COMPONENT_DOCUMENTATION)

  MESSAGE(STATUS "")
  MESSAGE(STATUS "Setting up documentation")

  ADD_SUBDIRECTORY(news)
  ADD_SUBDIRECTORY(doxygen)

  #
  # Install the static elements of the html documentation:
  #
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_DOCHTML_RELDIR}
    COMPONENT documentation
    PATTERN "CMakeLists.txt" EXCLUDE
    PATTERN "doxygen" EXCLUDE
    PATTERN "news" EXCLUDE
    )

  MESSAGE(STATUS "Setting up documentation - Done")
  MESSAGE(STATUS "")

ENDIF(DEAL_II_COMPONENT_DOCUMENTATION)

#
# Always install a minimalistic README and LICENSE file:
#

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/README.md
  ${CMAKE_SOURCE_DIR}/LICENSE.md
  DESTINATION ${DEAL_II_DOCREADME_RELDIR}
  COMPONENT library
  )

#
# Install summary.log an detailed.log
#

INSTALL(FILES
  ${CMAKE_BINARY_DIR}/summary.log
  ${CMAKE_BINARY_DIR}/detailed.log
  DESTINATION ${DEAL_II_DOCREADME_RELDIR}
  COMPONENT library
  )

#
# Add a dummy target to make documentation files known to IDEs.
#

FILE(GLOB _misc
  ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/headers/*.h
  )

ADD_LIBRARY(doxygen_headers OBJECT ${_misc})
SET_TARGET_PROPERTIES(doxygen_headers PROPERTIES LINKER_LANGUAGE C)
