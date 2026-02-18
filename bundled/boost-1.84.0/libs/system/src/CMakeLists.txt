## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2014 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(src_boost_system
  error_code.cpp
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
define_object_library(bundled_boost_system OBJECT ${src_boost_system})

