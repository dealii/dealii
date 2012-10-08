#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

FIND_PACKAGE(Perl REQUIRED)
FIND_PACKAGE(Doxygen REQUIRED)

IF(NOT DOXYGEN_FOUND OR NOT DOXYGEN_DOT_FOUND)
  # Just to be sure that we actually have dot...
  MESSAGE(FATAL_ERROR
    "Could not find doxygen and dot which is required for building the documentation"
    )
ENDIF()


#
# Prepare the example steps for inclusion:
#

ADD_SUBDIRECTORY(tutorial)


#
# Prepare auxiliary files for doxygen:
#

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/header.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/header.html
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/footer.html
  ${CMAKE_CURRENT_BINARY_DIR}/footer.html
  COPYONLY
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.css
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II/deal.css
  COPYONLY
  )

ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/images/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/images/structure.cmapx
  COMMAND ${DOXYGEN_DOT_EXECUTABLE}
  ARGS
    -Tpng -o ${CMAKE_CURRENT_BINARY_DIR}/images/structure.png
    -Tcmapx -o ${CMAKE_CURRENT_BINARY_DIR}/images/strucutre.cmapx
    ${CMAKE_CURRENT_SOURCE_DIR}/images/structure.dot
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/images/structure.dot
  )


#
# Finalize the doxygen configuration:
#

file(GLOB doxygen_input
  ${CMAKE_SOURCE_DIR}/include/deal.II/*
  ${CMAKE_BINARY_DIR}/include/deal.II/*
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/headers
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial
  )

TO_STRING(doxygen_input ${doxygen_input})

file(GLOB doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  ${CMAKE_SOURCE_DIR}/examples/*/doc
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/images
  )
TO_STRING(doxygen_image_path ${doxygen_image_path})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${doxygen_input}
  IMAGE_PATH=${doxygen_image_path}
  "
  )

ADD_CUSTOM_TARGET(doxygen ALL # TODO: Really all?
  #
  # ./deal.II is already created due to CONFIGURE_FILE(... deal.css ...)
  #
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/images/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/images/structure.cmapx
  COMMENT "Generating documentation via doxygen (This may take a _really_ long time)"
  VERBATIM
  )

ADD_DEPENDENCIES(doxygen tutorial)

INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}/doxygen
  COMPONENT documentation
  )

