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
# A glob. I'm sorry
#
file(GLOB DEAL_II_STEPS
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )

#
# Prepare the example steps for doxygen:
#
ADD_SUBDIRECTORY(tutorial) # uses ${DEAL_II_STEPS}
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
file(GLOB _doxygen_input
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/main.cpp
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/**/*.h

  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )
FOREACH(_step ${DEAL_II_STEPS})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)
  LIST(APPEND _doxygen_input
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
ENDFOREACH()
TO_STRING(_doxygen_input_string ${_doxygen_input})

file(GLOB _doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  ${CMAKE_SOURCE_DIR}/examples/*/doc
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/images

  ${CMAKE_CURRENT_BINARY_DIR}/images
  )
TO_STRING(_doxygen_image_path_string ${_doxygen_image_path})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${_doxygen_input_string}
  IMAGE_PATH=${_doxygen_image_path_string}
  "
  )

#
# And, finally, call doxygen:
#
ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_BINARY_DIR}/doxygen.log
  COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    > ${CMAKE_BINARY_DIR}/doxygen.log 2>&1 # *pssst*
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    tutorial
    images
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${_doxygen_input}
  COMMENT "Generating documentation via doxygen. This may take a _really_ long time..."
  VERBATIM
  )
ADD_CUSTOM_TARGET(doxygen ALL
  DEPENDS ${CMAKE_BINARY_DIR}/doxygen.log
  )


INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/deal.tag
  ${CMAKE_CURRENT_SOURCE_DIR}/stylesheet.css
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

