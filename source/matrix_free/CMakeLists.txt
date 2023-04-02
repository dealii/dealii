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


set(_inst
  evaluation_template_factory.inst.in
  evaluation_template_face_factory.inst.in
  evaluation_template_factory_hanging_nodes.inst.in
  fe_point_evaluation.inst.in
  mapping_info.inst.in
  matrix_free.inst.in
  shape_info.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/matrix_free/*.h
  )

if(DEAL_II_WITH_CUDA)
  set(_src
    cuda_matrix_free.cc
    ${_src}
    )
  set(_inst
    ${_inst}
    )
endif()

define_object_library(object_matrix_free OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_matrix_free "${_inst}")
