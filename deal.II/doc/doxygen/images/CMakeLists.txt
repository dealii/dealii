## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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
ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
  COMMAND ${DOXYGEN_DOT_EXECUTABLE}
    -Tpng -o ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    -Tcmapx -o ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
    ${CMAKE_CURRENT_SOURCE_DIR}/structure.dot
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/structure.dot
  )

ADD_CUSTOM_TARGET(images
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
  )
