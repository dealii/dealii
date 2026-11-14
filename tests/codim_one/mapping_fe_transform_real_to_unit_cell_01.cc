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


// Test MappingFE::transform_real_to_unit_cell() on curved cells of codimension
// one: map the quadrature points of the cell and of its faces to real space
// with transform_unit_to_real_cell() and back, which has to reproduce the
// points. The meshes are the surface of the octahedron with triangles and of
// the cube with quadrilaterals, both on the unit sphere, and the unit circle,
// with a spherical manifold on all cells.

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(const Triangulation<dim, spacedim> &tria,
     const FiniteElement<dim, spacedim> &fe_mapping)
{
  const MappingFE<dim, spacedim> mapping(fe_mapping);
  const ReferenceCell<dim>       reference_cell = fe_mapping.reference_cell();

  std::vector<Point<dim>> unit_points =
    reference_cell.get_gauss_type_quadrature(3).get_points();
  for (const unsigned int f : reference_cell.face_indices())
    {
      const Quadrature<dim> quadrature_on_face =
        QProjector<dim>::project_to_face(
          reference_cell,
          reference_cell.face_reference_cell(f).get_gauss_type_quadrature(3),
          f,
          numbers::default_geometric_orientation);
      unit_points.insert(unit_points.end(),
                         quadrature_on_face.get_points().begin(),
                         quadrature_on_face.get_points().end());
    }

  unsigned int n_points = 0, n_failed = 0;
  for (const auto &cell : tria.active_cell_iterators())
    for (const Point<dim> &p_unit : unit_points)
      {
        ++n_points;
        const Point<spacedim> p_real =
          mapping.transform_unit_to_real_cell(cell, p_unit);
        try
          {
            n_failed += mapping.transform_real_to_unit_cell(cell, p_real)
                          .distance(p_unit) > 1e-8;
          }
        catch (const typename Mapping<dim, spacedim>::ExcTransformationFailed &)
          {
            ++n_failed;
          }
      }

  deallog << fe_mapping.get_name() << ": failed round trips: " << n_failed
          << " of " << n_points << " points" << std::endl;
}



int
main()
{
  initlog();

  {
    // Octahedron with the outer normal as orientation of the triangles.
    const std::vector<Point<3>> vertices = {Point<3>(1, 0, 0),
                                            Point<3>(-1, 0, 0),
                                            Point<3>(0, 1, 0),
                                            Point<3>(0, -1, 0),
                                            Point<3>(0, 0, 1),
                                            Point<3>(0, 0, -1)};
    std::vector<CellData<2>>    cells;
    for (const unsigned int v0 : {0, 1})
      for (const unsigned int v1 : {2, 3})
        for (const unsigned int v2 : {4, 5})
          {
            CellData<2> cell(3);
            cell.vertices = {v0, v1, v2};
            if (cross_product_3d(vertices[v1] - vertices[v0],
                                 vertices[v2] - vertices[v0]) *
                  (vertices[v0] + vertices[v1] + vertices[v2]) <
                0)
              std::swap(cell.vertices[1], cell.vertices[2]);
            cells.push_back(cell);
          }

    Triangulation<2, 3> tria;
    tria.create_triangulation(vertices, cells, SubCellData());
    tria.set_all_manifold_ids(0);
    tria.set_manifold(0, SphericalManifold<2, 3>());
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_SimplexP<2, 3>(degree));
  }
  {
    Triangulation<2, 3> tria;
    GridGenerator::hyper_sphere(tria);
    tria.set_manifold(0, SphericalManifold<2, 3>());
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_Q<2, 3>(degree));
  }
  {
    Triangulation<1, 2> tria;
    GridGenerator::hyper_sphere(tria);
    tria.set_manifold(0, SphericalManifold<1, 2>());
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_Q<1, 2>(degree));
    for (unsigned int degree = 1; degree <= 3; ++degree)
      test(tria, FE_SimplexP<1, 2>(degree));
  }
}
