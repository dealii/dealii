## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(_src)
set(_inst)

if(DEAL_II_WITH_SYMENGINE)
  set(_src
    symengine_math.cc
    symengine_number_types.cc
    symengine_number_visitor_internal.cc
    symengine_optimizer.cc
    symengine_scalar_operations.cc
    symengine_tensor_operations.cc
    symengine_types.cc
    symengine_utilities.cc
    )

  set(_inst
    symengine_number_visitor_internal.inst.in
    symengine_optimizer.inst.in
    symengine_tensor_operations.inst.in
    )
endif()


file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/differentiation/sd/*.h
  )

define_object_library(object_differentiation_sd OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_differentiation_sd "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
