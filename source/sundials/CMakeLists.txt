## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


set(_src)
set(_inst)

if (DEAL_II_WITH_SUNDIALS)
  set(_src
    ${_src}
    arkode.cc
    ida.cc
    kinsol.cc
    n_vector.cc
    sunlinsol_wrapper.cc
    )

  set(_inst
    ${_inst}
    n_vector.inst.in
    )
endif()

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/sundials/*.h
  )

define_object_library(object_sundials OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_sundials "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
