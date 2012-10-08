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


#
# Prepare steps.png and steps.cmapx:
#

file(GLOB steps_args
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/steps.dot
  COMMAND ${PERL_EXECUTABLE}
  ARGS
    ${CMAKE_SOURCE_DIR}/doc/scripts/steps.pl
    ${steps_args}
    > ${CMAKE_CURRENT_BINARY_DIR}/steps.dot
  )

ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/steps.png
    ${CMAKE_CURRENT_BINARY_DIR}/steps.cmapx
  COMMAND ${DOXYGEN_DOT_EXECUTABLE}
  ARGS
    -Tpng -o ${CMAKE_CURRENT_BINARY_DIR}/steps.png
    -Tcmapx -o ${CMAKE_CURRENT_BINARY_DIR}/steps.cmapx
    ${CMAKE_CURRENT_BINARY_DIR}/steps.dot
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/steps.dot
  )


#
# Prepare toc.html:
#

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/toc.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/toc.html.in
  COPYONLY
  )

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/toc.html
  COMMAND ${PERL_EXECUTABLE}
  ARGS
    ${CMAKE_SOURCE_DIR}/doc/scripts/make_toc.pl
    > ${CMAKE_CURRENT_BINARY_DIR}/toc.html
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/toc.html.in
    ${CMAKE_CURRENT_BINARY_DIR}/steps.cmapx
  )

ADD_CUSTOM_TARGET(blubb ALL
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/toc.html
  )


#
# Prepare toc.html:
#

