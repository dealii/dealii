## ---------------------------------------------------------------------
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
# Generate the source documentation via doxygen:
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
IF(NOT DOXYGEN_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find doxygen which is required for building the documentation"
    )
ENDIF()

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

#
# Prepare auxiliary files for doxygen:
#
IF(DEAL_II_DOXYGEN_USE_MATHJAX)
  SET(_use_mathjax YES)
ELSE()
  SET(_use_mathjax NO)
ENDIF()
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mod_header.pl.in
  ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mod_footer.pl.in
  ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl
  )


#
# Generate header, footer and style files for doxygen
#
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/header.html
    ${CMAKE_CURRENT_BINARY_DIR}/footer.html
  COMMAND ${DOXYGEN_EXECUTABLE} -w html header.html footer.html style.css options.dox
  COMMAND ${PERL_EXECUTABLE} -pi~ ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl header.html
  COMMAND ${PERL_EXECUTABLE} -pi~ ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl footer.html
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl
  )

#
# Finalize the doxygen configuration:
#

SET(_doxygen_input
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/
)

LIST(APPEND _doxygen_input
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/doc/news/
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/
  )

SET(_doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/images
  ${CMAKE_CURRENT_BINARY_DIR}/images
  )

file(GLOB _doxygen_depend
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/main.cpp
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/**/*.h
  )
LIST(APPEND _doxygen_depend
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )
FOREACH(_step ${DEAL_II_STEPS})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)
  LIST(APPEND _doxygen_depend
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
  LIST(APPEND _doxygen_input
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
ENDFOREACH()

TO_STRING(_doxygen_image_path_string ${_doxygen_image_path})
TO_STRING(_doxygen_input_string ${_doxygen_input})

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
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/header.html
    ${CMAKE_CURRENT_BINARY_DIR}/footer.html
    ${_doxygen_depend}
  COMMENT "Generating documentation via doxygen."
  VERBATIM
  )
ADD_CUSTOM_TARGET(doxygen ALL
  DEPENDS ${CMAKE_BINARY_DIR}/doxygen.log
  )
ADD_DEPENDENCIES(documentation doxygen)

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/deal.tag
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/doc/deal.ico
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen/deal.II
  COMPONENT documentation
  )

INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

