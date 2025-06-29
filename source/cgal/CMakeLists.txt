## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2022 - 2025 by the deal.II authors
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

if(DEAL_II_WITH_CGAL)
  #
  # We have to compile the "intersections.cc" file without the misleading
  # indentation warning enabled. Otherwise, we run into quite a number of
  # warnings with gcc, see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89549
  #
  enable_if_supported(_flag "-Wno-misleading-indentation")
  set_property(SOURCE "intersections.cc"
    APPEND PROPERTY COMPILE_OPTIONS "${_flag}"
    )

  set(_src
    ${_src}
    surface_mesh.cc
    intersections.cc
    polygon.cc
  )

  set(_inst
      ${_inst}
  surface_mesh.inst.in
  intersections.inst.in
  polygon.inst.in
  )
endif()

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/cgal/*.h
  )

define_object_library(object_cgal OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_cgal "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
