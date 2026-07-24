// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check TriaAccessor::real_to_unit_cell_affine_approximation for non-hypercube
// reference cells.

#include <deal.II/base/derivative_form.h>

#include <deal.II/fe/mapping_p1.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include <string>
#include <vector>


template <int dim, int spacedim>
DerivativeForm<1, dim, spacedim>
make_linear_part()
{
  DerivativeForm<1, dim, spacedim> A;
  for (unsigned int d = 0; d < spacedim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      A[d][e] = (d == e ? 1.0 + 0.2 * d : 0.0) +
                0.05 * (d + 1) * (e + 2);

  return A;
}



template <int spacedim>
Tensor<1, spacedim>
make_offset()
{
  Tensor<1, spacedim> b;
  for (unsigned int d = 0; d < spacedim; ++d)
    b[d] = -0.3 + 0.2 * d;

  return b;
}



template <int dim, int spacedim>
Point<spacedim>
transform(const DerivativeForm<1, dim, spacedim> &A,
          const Tensor<1, spacedim>              &b,
          const Point<dim>                       &unit_point)
{
  Point<spacedim> real_point;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      real_point[d] = b[d];
      for (unsigned int e = 0; e < dim; ++e)
        real_point[d] += A[d][e] * unit_point[e];
    }

  return real_point;
}



template <int dim, int spacedim>
void
test(const ReferenceCell<dim>      &reference_cell,
     const std::vector<Point<dim>> &unit_points,
     const std::string             &name,
     const bool                     compare_with_mapping_p1)
{
  const auto A = make_linear_part<dim, spacedim>();
  const auto b = make_offset<spacedim>();

  std::vector<Point<spacedim>> vertices(reference_cell.n_vertices());
  for (unsigned int v = 0; v < reference_cell.n_vertices(); ++v)
    vertices[v] = transform(A, b, reference_cell.vertex(v));

  std::vector<CellData<dim>> cells(
    1, CellData<dim>(reference_cell.n_vertices()));
  for (unsigned int v = 0; v < reference_cell.n_vertices(); ++v)
    cells[0].vertices[v] = v;

  Triangulation<dim, spacedim> triangulation;
  triangulation.create_triangulation(vertices, cells, SubCellData());

  const auto cell = triangulation.begin_active();

  MappingP1<dim, spacedim> mapping_p1;
  for (const Point<dim> &unit_point : unit_points)
    {
      const Point<spacedim> real_point = transform(A, b, unit_point);
      const Point<dim>      approximate =
        cell->real_to_unit_cell_affine_approximation(real_point);

      AssertThrow(approximate.distance(unit_point) < 1e-11,
                  ExcMessage("Unexpected real-to-unit affine approximation."));

      if (compare_with_mapping_p1)
        {
          const Point<dim> mapping_p1_point =
            mapping_p1.transform_real_to_unit_cell(cell, real_point);
          AssertThrow(approximate.distance(mapping_p1_point) < 1e-11,
                      ExcMessage("Unexpected difference from MappingP1."));
        }
    }

  deallog << name << " OK" << std::endl;
}



int
main()
{
  initlog();

  const std::vector<Point<2>> triangle_points = {
    ReferenceCells::Triangle.vertex(0),
    ReferenceCells::Triangle.barycenter(),
    Point<2>{0.25, 0.25},
    Point<2>{0.6, 0.2}};
  const std::vector<Point<3>> tetrahedron_points = {
    ReferenceCells::Tetrahedron.vertex(0),
    ReferenceCells::Tetrahedron.barycenter(),
    Point<3>{0.1, 0.2, 0.3},
    Point<3>{0.3, 0.2, 0.1}};
  const std::vector<Point<3>> wedge_points = {
    ReferenceCells::Wedge.vertex(0),
    ReferenceCells::Wedge.barycenter(),
    Point<3>{0.2, 0.3, 0.4},
    Point<3>{0.1, 0.1, 0.8}};
  const std::vector<Point<3>> pyramid_points = {
    ReferenceCells::Pyramid.vertex(0),
    ReferenceCells::Pyramid.barycenter(),
    Point<3>{0.1, -0.1, 0.5},
    Point<3>{-0.2, 0.2, 0.3}};

  test<2, 2>(ReferenceCells::Triangle,
             triangle_points,
             "Triangle",
             true);
  test<3, 3>(ReferenceCells::Tetrahedron,
             tetrahedron_points,
             "Tetrahedron",
             true);
  test<2, 3>(ReferenceCells::Triangle,
             triangle_points,
             "Triangle codimension one",
             true);
  test<3, 3>(ReferenceCells::Wedge,
             wedge_points,
             "Wedge",
             false);
  test<3, 3>(ReferenceCells::Pyramid,
             pyramid_points,
             "Pyramid",
             false);
}
