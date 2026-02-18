## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(_src)
set(_inst)

if(DEAL_II_TRILINOS_WITH_NOX)
  set(_src
    ${_src}
    nox.cc
    )

  set(_inst
    ${_inst}
    nox.inst.in
    )
endif()


file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/trilinos/*.h
  )

define_object_library(object_trilinos OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_trilinos "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
