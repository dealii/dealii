## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# All source/base files are sufficiently small that they may be added to the
# unity build file (see the developer documentation on DEAL_II_UNITY_BUILD
# for more information).
#
set(_unity_include_src
  all.cc
  )

# Add CUDA wrapper files
if(DEAL_II_WITH_CUDA)
  set(_separate_src
    ${_separate_src}
    cuda.cc
  )
endif()

# determined by profiling
set(_n_includes_per_unity_file 29)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  bounding_box.inst.in
  data_out_base.inst.in
  function.inst.in
  function_signed_distance.inst.in
  function_restriction.inst.in
  function_time.inst.in
  function_tools.inst.in
  incremental_function.inst.in
  geometric_utilities.inst.in
  mpi.inst.in
  mpi_noncontiguous_partitioner.inst.in
  mpi_remote_point_evaluation.inst.in
  partitioner.inst.in
  polynomials_rannacher_turek.inst.in
  symmetric_tensor.inst.in
  symbolic_function.inst.in
  tensor_function.inst.in
  tensor_function_parser.inst.in
  mu_parser_internal.inst.in
  time_stepping.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.h
  )

define_object_library(object_base OBJECT ${_src} ${_header} ${_inst}
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )
expand_instantiations(object_base "${_inst}")
