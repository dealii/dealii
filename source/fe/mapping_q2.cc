// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q2.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingQ2<dim, spacedim>::MappingQ2(
  std::vector<std::vector<Point<spacedim>>> &support_points)
  : MappingQGeneric<dim, spacedim>(2)
  , support_points(support_points)
{}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ2<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<MappingQ2<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::vector<Point<spacedim>>
MappingQ2<dim, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  // check if level is 0, since currently, no mesh refinement is supported
  AssertDimension(cell->level(), 0);

  // number of total points is 9 (3^2) for quad9 element
  // and 27 (3^3) for hex27 element
  const unsigned int num_points = Utilities::pow(3, dim);

  // vector that will hold the vertices plus support points
  std::vector<Point<spacedim>> result(num_points);

  // get vertices of the current cell and add to results
  const auto &vertices = this->get_vertices(cell);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    result[i] = vertices[i];

  // get vector with support points for current cell and add to results
  const std::vector<Point<spacedim>> &sup_points =
    support_points[cell->active_cell_index()];
  for (unsigned int i = 0;
       i < num_points - GeometryInfo<dim>::vertices_per_cell;
       ++i)
    result[i + GeometryInfo<dim>::vertices_per_cell] = sup_points[i];

  return result;

  // new version since now vertices are also stored in the support points vector
  //  return support_points[cell->active_cell_index()];
}



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q2.inst"



DEAL_II_NAMESPACE_CLOSE
