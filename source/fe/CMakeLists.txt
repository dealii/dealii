## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  block_mask.cc
  component_mask.cc
  fe_abf.cc
  fe_bdm.cc
  fe.cc
  fe_data.cc
  fe_dgp.cc
  fe_dgp_monomial.cc
  fe_dgp_nonparametric.cc
  fe_dgq.cc
  fe_dg_vector.cc
  fe_face.cc
  fe_nedelec.cc
  fe_nothing.cc
  fe_poly.cc
  fe_poly_tensor.cc
  fe_q_base.cc
  fe_q.cc
  fe_q_dg0.cc
  fe_q_hierarchical.cc
  fe_q_iso_q1.cc
  fe_raviart_thomas.cc
  fe_raviart_thomas_nodal.cc
  fe_system.cc
  fe_tools.cc
  fe_tools_interpolate.cc
  fe_trace.cc
  fe_values.cc
  fe_values_inst2.cc
  mapping_c1.cc
  mapping_cartesian.cc
  mapping.cc
  mapping_q1.cc
  mapping_q1_eulerian.cc
  mapping_q.cc
  mapping_q_eulerian.cc
  )

SET(_inst
  fe_abf.inst.in
  fe_bdm.inst.in
  fe_dgp.inst.in
  fe_dgp_monomial.inst.in
  fe_dgp_nonparametric.inst.in
  fe_dgq.inst.in
  fe_dg_vector.inst.in
  fe_face.inst.in
  fe.inst.in
  fe_nedelec.inst.in
  fe_nothing.inst.in
  fe_poly.inst.in
  fe_poly_tensor.inst.in
  fe_q_base.inst.in
  fe_q_dg0.inst.in
  fe_q_hierarchical.inst.in
  fe_q.inst.in
  fe_q_iso_q1.inst.in
  fe_raviart_thomas.inst.in
  fe_raviart_thomas_nodal.inst.in
  fe_system.inst.in
  fe_tools.inst.in
  fe_tools_interpolate.inst.in
  fe_trace.inst.in
  fe_values.decl.1.inst.in
  fe_values.decl.2.inst.in
  fe_values.impl.1.inst.in
  fe_values.impl.2.inst.in
  fe_values.inst.in
  mapping_c1.inst.in
  mapping_cartesian.inst.in
  mapping.inst.in
  mapping_q1_eulerian.inst.in
  mapping_q1.inst.in
  mapping_q_eulerian.inst.in
  mapping_q.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/fe/*.h
  )

DEAL_II_ADD_LIBRARY(obj_fe OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_fe "${_inst}")
