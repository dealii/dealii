## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(_src)
set(_inst)

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/optimization/*.h
  )

if(DEAL_II_WITH_PETSC)
  list(APPEND _src
    petsc_linear_optimizer.cc
    petsc_nonlinear_optimizer.cc
    )
endif()

define_object_library(object_optimization OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_optimization "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")

add_subdirectory(rol)
