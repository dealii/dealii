## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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

IF(DEAL_II_COMPONENT_PYTHON_BINDINGS)

  MESSAGE(STATUS "")
  MESSAGE(STATUS "Setting up python bindings")

  ADD_SUBDIRECTORY(source)

  MESSAGE(STATUS "Setting up python bindings - Done")
  MESSAGE(STATUS "")

ENDIF(DEAL_II_COMPONENT_PYTHON_BINDINGS)
