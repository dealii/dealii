## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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

########################################################################
#
# Process the tutorial and code-gallery files into inputs for doxygen
#
########################################################################

ADD_SUBDIRECTORY(tutorial)
ADD_SUBDIRECTORY(code-gallery)


########################################################################
#
# Set up all of the other input pieces we want to give to doxygen
#
########################################################################

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
  ${CMAKE_SOURCE_DIR}/source/
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/doc/news/
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# Add other directories in which to find images
SET(_doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  )

file(GLOB _doxygen_depend
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/**/*.h
  )
LIST(APPEND _doxygen_depend
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# find all tutorial programs so we can add dependencies as appropriate
FILE(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )
FOREACH(_step ${_deal_ii_steps})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)
  LIST(APPEND _doxygen_depend
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
  LIST(APPEND _doxygen_input
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
ENDFOREACH()

# Also find all code gallery programs (if available) for the same reason.
# The logic here follows the same as in code-gallery/CMakeLists.txt
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
IF (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)
  FILE(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  FOREACH(_step ${_code_gallery_names})
    GET_FILENAME_COMPONENT(_step "${_step}" NAME)

    LIST(APPEND _doxygen_depend
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
    LIST(APPEND _doxygen_input
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
  ENDFOREACH()
ENDIF()


# always make the doxygen run depend on the code-gallery.h file
# (whether generated from the code gallery or copied from
# no-code-gallery.h; both happen in code-gallery/CMakeLists.txt)
LIST(APPEND _doxygen_input
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)
LIST(APPEND _doxygen_depend
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)


TO_STRING(_doxygen_image_path_string ${_doxygen_image_path})
TO_STRING(_doxygen_input_string ${_doxygen_input})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${_doxygen_input_string}
  IMAGE_PATH=${_doxygen_image_path_string}
  "
  )

########################################################################
#
# And, finally, call doxygen:
#
########################################################################

ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_BINARY_DIR}/doxygen.log
  COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    > ${CMAKE_BINARY_DIR}/doxygen.log 2>&1 # *pssst*
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    tutorial
    code-gallery
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/header.html
    ${CMAKE_CURRENT_BINARY_DIR}/footer.html
    ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter
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

