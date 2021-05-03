## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2014 - 2019 by the deal.II authors
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

INCLUDE_DIRECTORIES(
  ${BOOST_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

SET(src_boost_system
  error_code.cpp
  )

DEAL_II_ADD_LIBRARY(obj_boost_system OBJECT ${src_boost_system})

