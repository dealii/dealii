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


// Test MappingFE::transform_real_to_unit_cell() on curved simplex cells: map
// points of the reference cell to real space with
// transform_unit_to_real_cell() and back, which has to reproduce the points.
// The meshes are converted from hypercube meshes with a curved manifold on the
// boundary only, i.e., the cells at the boundary have curved faces and
// straight interior edges. A curved manifold on all cells would invert some
// of the cells at the center of the ball and at the axis of the cylinder.

#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
convert(Triangulation<dim>  &tria_hypercube,
        Triangulation<dim>  &tria,
        const Manifold<dim> &manifold)
{
  tria_hypercube.refine_global(1);
  Triangulation<dim> tria_flat;
  GridGenerator::flatten_triangulation(tria_hypercube, tria_flat);
  tria_flat.set_all_manifold_ids(numbers::flat_manifold_id);
  tria_flat.set_all_manifold_ids_on_boundary(0 /* boundary id */, 0);
  tria_flat.set_manifold(0, manifold);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_flat, tria);
  tria.set_manifold(0, manifold);
}



template <int dim>
void
test(const std::string &name, Triangulation<dim> &tria)
{
  const ReferenceCell<dim> reference_cell = ReferenceCells::get_simplex<dim>();
  const QGaussSimplex<dim - 1> quadrature_face(3);
  const QGaussSimplex<dim>     quadrature_cell(3);

  for (unsigned int mapping_degree = 2; mapping_degree <= 3; ++mapping_degree)
    {
      const MappingFE<dim> mapping(FE_SimplexP<dim>{mapping_degree});

      const FE_Nothing<dim> fe(reference_cell);
      FEValues<dim>         fe_values(mapping,
                              fe,
                              QWitherdenVincentSimplex<dim>(5),
                              update_jacobians);

      unsigned int n_face = 0, n_failed_face = 0, n_cell = 0, n_failed_cell = 0;

      const auto failed = [&](const auto &cell, const Point<dim> &p_unit) {
        const Point<dim> p_real =
          mapping.transform_unit_to_real_cell(cell, p_unit);
        try
          {
            return mapping.transform_real_to_unit_cell(cell, p_real)
                     .distance(p_unit) > 1e-8;
          }
        catch (const typename Mapping<dim>::ExcTransformationFailed &)
          {
            return true;
          }
      };

      for (const auto &cell : tria.active_cell_iterators())
        if (cell->at_boundary())
          {
            // All cells have to be valid, i.e., the map is invertible.
            fe_values.reinit(cell);
            for (const unsigned int q : fe_values.quadrature_point_indices())
              AssertThrow(fe_values.jacobian(q).determinant() > 0,
                          ExcInternalError());

            for (const unsigned int f : cell->face_indices())
              if (cell->at_boundary(f))
                {
                  const Quadrature<dim> quadrature_on_face =
                    QProjector<dim>::project_to_face(
                      reference_cell,
                      quadrature_face,
                      f,
                      numbers::default_geometric_orientation);
                  for (const Point<dim> &p_unit :
                       quadrature_on_face.get_points())
                    {
                      ++n_face;
                      n_failed_face += failed(cell, p_unit);
                    }
                }
            for (const Point<dim> &p_unit : quadrature_cell.get_points())
              {
                ++n_cell;
                n_failed_cell += failed(cell, p_unit);
              }
          }

      deallog << name << ", mapping degree " << mapping_degree
              << ": failed round trips: " << n_failed_face << " of " << n_face
              << " face quadrature points on the boundary, " << n_failed_cell
              << " of " << n_cell << " cell quadrature points" << std::endl;
    }
}



int
main()
{
  initlog();

  {
    Triangulation<2> tria_hypercube, tria;
    GridGenerator::hyper_ball(tria_hypercube);
    convert(tria_hypercube, tria, SphericalManifold<2>());
    test("2D ball", tria);
  }
  {
    Triangulation<3> tria_hypercube, tria;
    GridGenerator::hyper_ball(tria_hypercube);
    convert(tria_hypercube, tria, SphericalManifold<3>());
    test("3D ball", tria);
  }
  {
    Triangulation<3> tria_hypercube, tria;
    GridGenerator::cylinder(tria_hypercube, 1.0, 1.0);
    convert(tria_hypercube, tria, CylindricalManifold<3>(0));
    test("3D cylinder", tria);
  }
}
