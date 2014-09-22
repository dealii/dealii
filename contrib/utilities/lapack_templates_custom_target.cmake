## ---------------------------------------------------------------------
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
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/
  COMMAND perl
  ARGS ${CMAKE_SOURCE_DIR}/scripts/lapack_templates.pl
       ${CMAKE_CURRENT_SOURCE_DIR}/deal.II/lac/lapack_templates.h.in
       > ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  )

ADD_CUSTOM_TARGET(lapack_templates ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/deal.II/lac/lapack_templates.h
  )
