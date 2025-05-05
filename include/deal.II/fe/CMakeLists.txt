## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------



#
# Keep track of all of the header files in this directory
#
set(_header_files
    block_mask.h
    component_mask.h
    fe_abf.h
    fe_base.h
    fe_bdm.h
    fe_bernardi_raugel.h
    fe_bernstein.h
    fe_coupling_values.h
    fe_data.h
    fe_dgp.h
    fe_dgp_monomial.h
    fe_dgp_nonparametric.h
    fe_dgq.h
    fe_dg_vector.h
    fe_dg_vector.templates.h
    fe_enriched.h
    fe_face.h
    fe.h
    fe_hermite.h
    fe_interface_values.h
    fe_nedelec.h
    fe_nedelec_sz.h
    fe_nothing.h
    fe_p1nc.h
    fe_poly_face.h
    fe_poly_face.templates.h
    fe_poly.h
    fe_poly_tensor.h
    fe_pyramid_p.h
    fe_q_base.h
    fe_q_bubbles.h
    fe_q_dg0.h
    fe_q.h
    fe_q_hierarchical.h
    fe_q_iso_q1.h
    fe_rannacher_turek.h
    fe_raviart_thomas.h
    fe_rt_bubbles.h
    fe_series.h
    fe_simplex_p_bubbles.h
    fe_simplex_p.h
    fe_system.h
    fe_tools_extrapolate.templates.h
    fe_tools.h
    fe_tools_interpolate.templates.h
    fe_tools.templates.h
    fe_trace.h
    fe_update_flags.h
    fe_values_base.h
    fe_values_extractors.h
    fe_values.h
    fe_values_views.h
    fe_values_views_internal.h
    fe_wedge_p.h
    mapping_c1.h
    mapping_cartesian.h
    mapping_fe_field.h
    mapping_fe.h
    mapping.h
    mapping_manifold.h
    mapping_q1_eulerian.h
    mapping_q1.h
    mapping_q_cache.h
    mapping_q_eulerian.h
    mapping_q_generic.h
    mapping_q.h
    mapping_q_internal.h
    mapping_related_data.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
