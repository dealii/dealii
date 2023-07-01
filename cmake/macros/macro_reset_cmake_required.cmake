## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2022 by the deal.II authors
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

#
# A small macro to reset the CMAKE_REQUIRED_* variables to its default
# values
#
# Usage:
#     RESET_CMAKE_REQUIRED_FLAGS
#

macro(reset_cmake_required)
  set(CMAKE_REQUIRED_FLAGS ${DEAL_II_CXX_FLAGS_SAVED})
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_LIBRARIES ${DEAL_II_LINKER_FLAGS_SAVED})
endmacro()

