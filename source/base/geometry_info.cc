// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
constexpr unsigned int GeometryInfo<dim>::max_children_per_cell;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::faces_per_cell;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::max_children_per_face;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::vertices_per_cell;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::vertices_per_face;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::lines_per_face;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::quads_per_face;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::lines_per_cell;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::quads_per_cell;
template <int dim>
constexpr unsigned int GeometryInfo<dim>::hexes_per_cell;

template <int dim>
constexpr std::array<int, GeometryInfo<dim>::faces_per_cell>
  GeometryInfo<dim>::unit_normal_orientation;

template <int dim>
constexpr std::array<std::array<unsigned int, dim>,
                     GeometryInfo<dim>::vertices_per_cell>
  GeometryInfo<dim>::vertex_to_face;

template <int dim>
constexpr std::array<unsigned int, GeometryInfo<dim>::faces_per_cell>
  GeometryInfo<dim>::unit_normal_direction;

template <int dim>
constexpr std::array<Tensor<1, dim>, GeometryInfo<dim>::faces_per_cell>
  GeometryInfo<dim>::unit_normal_vector;

template <int dim>
constexpr std::array<std::array<Tensor<1, dim>, dim - 1>,
                     GeometryInfo<dim>::faces_per_cell>
  GeometryInfo<dim>::unit_tangential_vectors;

template <int dim>
constexpr std::array<unsigned int, GeometryInfo<dim>::vertices_per_cell>
  GeometryInfo<dim>::dx_to_deal;

template <int dim>
constexpr std::array<unsigned int, GeometryInfo<dim>::faces_per_cell>
  GeometryInfo<dim>::opposite_face;

template <int dim>
constexpr std::array<unsigned int, GeometryInfo<dim>::vertices_per_cell>
  GeometryInfo<dim>::ucd_to_deal;

const std::array<unsigned int, GeometryInfo<0>::vertices_per_cell>
  GeometryInfo<0>::ucd_to_deal = {{0}};

const std::array<unsigned int, GeometryInfo<0>::vertices_per_cell>
  GeometryInfo<0>::dx_to_deal = {{0}};

template struct GeometryInfo<1>;
template struct GeometryInfo<2>;
template struct GeometryInfo<3>;
template struct GeometryInfo<4>;

template void
GeometryInfo<1>::alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<1> (&)[vertices_per_cell],
   Tensor<1 - 1, 1> (&)[vertices_per_cell])
#else
  (const Point<1> *, Tensor<1 - 1, 1> *)
#endif
    ;

template void
GeometryInfo<1>::alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<2> (&)[vertices_per_cell],
   Tensor<2 - 1, 2> (&)[vertices_per_cell])
#else
  (const Point<2> *, Tensor<2 - 1, 2> *)
#endif
    ;

template void
GeometryInfo<2>::alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<2> (&vertices)[vertices_per_cell],
   Tensor<2 - 2, 2> (&forms)[vertices_per_cell])
#else
  (const Point<2> *, Tensor<2 - 2, 2> *)
#endif
    ;

template void
GeometryInfo<2>::alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<3> (&vertices)[vertices_per_cell],
   Tensor<3 - 2, 3> (&forms)[vertices_per_cell])
#else
  (const Point<3> *, Tensor<3 - 2, 3> *)
#endif
    ;


template void
GeometryInfo<3>::alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<3> (&vertices)[vertices_per_cell],
   Tensor<3 - 3, 3> (&forms)[vertices_per_cell])
#else
  (const Point<3> *, Tensor<3 - 3, 3> *)
#endif
    ;

DEAL_II_NAMESPACE_CLOSE
