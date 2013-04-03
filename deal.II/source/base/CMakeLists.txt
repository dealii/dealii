#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.h
  )

SET(_src
  auto_derivative_function.cc
  conditional_ostream.cc
  config.cc
  convergence_table.cc
  data_out_base.cc
  # do NOT include dummy.cc
  event.cc
  exceptions.cc
  flow_function.cc
  function.cc
  function_derivative.cc
  function_lib.cc
  function_lib_cutoff.cc
  function_parser.cc
  function_time.cc
  geometry_info.cc
  index_set.cc
  job_identifier.cc
  logstream.cc
  memory_consumption.cc
  mpi.cc
  multithread_info.cc
  parallel.cc
  parameter_handler.cc
  parsed_function.cc
  partitioner.cc
  path_search.cc
  polynomial.cc
  polynomials_abf.cc
  polynomials_adini.cc
  polynomials_bdm.cc
  polynomials_nedelec.cc
  polynomial_space.cc
  polynomials_p.cc
  polynomials_raviart_thomas.cc
  quadrature.cc
  quadrature_lib.cc
  quadrature_selector.cc
  subscriptor.cc
  symmetric_tensor.cc
  table_handler.cc
  tensor_function.cc
  tensor_product_polynomials.cc
  tensor_product_polynomials_const.cc
  thread_management.cc
  timer.cc
  utilities.cc
  )

DEAL_II_ADD_LIBRARY(obj_base OBJECT ${_src} ${_header}
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )

