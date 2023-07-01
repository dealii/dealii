## ---------------------------------------------------------------------
##
## Copyright (C) 2022 - 2023 by the deal.II authors
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
# We have to compile the "intersections.cc" file without the misleading
# indentation warning enabled. Otherwise, we run into quite a number of
# warnings with gcc, see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89549
#
if(DEAL_II_WITH_CGAL)
  enable_if_supported(_flag "-Wno-misleading-indentation")
  set_property(SOURCE "intersections.cc"
    APPEND PROPERTY COMPILE_OPTIONS "${_flag}"
    )
endif()

set(_src
  surface_mesh.cc
  intersections.cc
)

set(_inst
surface_mesh.inst.in
intersections.inst.in
)

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/cgal/*.h
  )

define_object_library(object_cgal OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_cgal "${_inst}")
