## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

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
