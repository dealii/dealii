## ------------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Set up all necessary bits for the documentation
#
if(DEAL_II_COMPONENT_DOCUMENTATION)

  message(STATUS "")
  message(STATUS "Setting up documentation")

  add_subdirectory(news)
  add_subdirectory(doxygen)

  #
  # Install the static elements of the html documentation:
  #
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_DOCHTML_RELDIR}
    COMPONENT documentation
    PATTERN "CMakeLists.txt" EXCLUDE
    PATTERN "doxygen" EXCLUDE
    PATTERN "news" EXCLUDE
    )

  message(STATUS "Setting up documentation - Done")
  message(STATUS "")

endif()

#
# Always install a minimalistic README and LICENSE file:
#

install(FILES
  ${CMAKE_SOURCE_DIR}/README.md
  ${CMAKE_SOURCE_DIR}/LICENSE.md
  DESTINATION ${DEAL_II_DOCREADME_RELDIR}
  COMPONENT library
  )

#
# Install summary.log an detailed.log
#

install(FILES
  ${CMAKE_BINARY_DIR}/summary.log
  ${CMAKE_BINARY_DIR}/detailed.log
  DESTINATION ${DEAL_II_DOCREADME_RELDIR}
  COMPONENT library
  )

#
# Add a dummy target to make documentation files known to IDEs.
#

file(GLOB _misc
  ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/headers/*.h
  )

add_library(doxygen_headers OBJECT ${_misc})
set_target_properties(doxygen_headers PROPERTIES LINKER_LANGUAGE C)
