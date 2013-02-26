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
# Prepare the example steps for doxygen:
#

ADD_SUBDIRECTORY(tutorial)
ADD_SUBDIRECTORY(images)


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

#
# Finalize the doxygen configuration:
#

file(GLOB doxygen_input
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/main.cpp
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/**/*.h

  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial
  )
TO_STRING(doxygen_input_string ${doxygen_input})

file(GLOB doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/images

  ${CMAKE_CURRENT_BINARY_DIR}/images
  )
TO_STRING(doxygen_image_path_string ${doxygen_image_path})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${doxygen_input_string}
  IMAGE_PATH=${doxygen_image_path_string}
  "
  )

#
# And, finally, call doxygen:
#

ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.log
  COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    > ${CMAKE_BINARY_DIR}/doxygen.log 2>&1 # *pssst*
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    tutorial
    images
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${doxygen_input}
  COMMENT "Generating documentation via doxygen (This may take a _really_ long time)"
  VERBATIM
  )
ADD_CUSTOM_TARGET(doxygen ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.log
  )

INSTALL(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/stylesheet.css
  DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}/doxygen
  COMPONENT documentation
  )

INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}/doxygen
  COMPONENT documentation
  )

