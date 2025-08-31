// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/non_matching/closest_surface_point.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/vector_tools.h>

#include <cmath>

#include "../tests.h"



// Test for NonMatching::ClosestSurfacePoint using a unit sphere level set


template <int dim>
void
test()
{
  deallog << "Testing ClosestSurfacePoint in " << dim << "D" << std::endl;

  // Create hypercube [-1.2, 1.2]^dim
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1.2, 1.2);

  // Refine 3 times
  triangulation.refine_global(3);

  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;

  // Create DoFHandler with FE_Q(4)
  DoFHandler<dim> dof_handler(triangulation);
  // Use high order FE so that the points are close to unit sphere
  FE_Q<dim> fe(4);
  dof_handler.distribute_dofs(fe);

  deallog << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  // Create level set vector and interpolate unit sphere function
  Vector<double>                               level_set(dof_handler.n_dofs());
  const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;

  // Interpolate the signed distance function
  VectorTools::interpolate(dof_handler, signed_distance_sphere, level_set);

  // Create ClosestSurfacePoint object
  typename NonMatching::ClosestSurfacePoint<dim, double>::AdditionalData data;
  data.tolerance    = 1e-10;
  data.n_iterations = 20;

  MappingCartesian<dim>                         mapping;
  NonMatching::ClosestSurfacePoint<dim, double> closest_point_finder(
    level_set, dof_handler, mapping, data);

  NonMatching::MeshClassifier<dim> mesh_classifier(dof_handler, level_set);
  mesh_classifier.reclassify();

  // Test parameters
  const double tolerance =
    dof_handler.begin_active()->minimum_vertex_distance() * 1e-2;
  unsigned int n_tested_points    = 0;
  unsigned int n_points_on_sphere = 0;

  // Limit the number of test cells to cut down on test time
  int n_test_cells = 20;

  // Loop over all cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // Check if cell potentially intersects the unit sphere
      // (cell contains points both inside and outside unit sphere)
      std::vector<Point<dim>> vertices;
      bool                    has_inside  = false;
      bool                    has_outside = false;

      const NonMatching::LocationToLevelSet cell_location =
        mesh_classifier.location_to_level_set(cell);
      if (cell_location == NonMatching::LocationToLevelSet::outside)
        continue;

      for (const auto &face : GeometryInfo<dim>::face_indices())
        {
          if (n_test_cells <= 0)
            continue;

          // Check if cell has a neighbor on this face
          if (cell->at_boundary(face))
            continue;

          const auto                            neighbor = cell->neighbor(face);
          const NonMatching::LocationToLevelSet neighbor_location =
            mesh_classifier.location_to_level_set(neighbor);

          // ignore neighbors that are inside the level set
          if (neighbor_location == NonMatching::LocationToLevelSet::inside)
            continue;

          // Decrease the number of remaining test cells
          --n_test_cells;

          // neighbor is outside or intersected, so we can test the face
          // Create test points on the face
          QGauss<dim - 1>   face_quadrature(3);
          FEFaceValues<dim> fe_face_values(fe,
                                           face_quadrature,
                                           update_quadrature_points);
          fe_face_values.reinit(cell, face);

          // Get the quadrature points on the face
          std::vector<Point<dim>> face_points;
          for (unsigned int q = 0; q < face_quadrature.size(); ++q)
            {
              face_points.push_back(fe_face_values.quadrature_point(q));
            }

          // Run closest point finder: search_cell=neighbor,
          // reference_cell=cell
          auto closest_points_result =
            closest_point_finder.compute_closest_surface_points(neighbor,
                                                                cell,
                                                                face_points);

          // Extract the results
          const auto &closest_real_points = closest_points_result.first;
          const auto &closest_unit_reference_points =
            closest_points_result.second;

          // Test the results
          for (unsigned int q = 0; q < face_points.size(); ++q)
            {
              n_tested_points++;

              const Point<dim> &original_point = face_points[q];
              const Point<dim> &closest_point  = closest_real_points[q];

              // Calculate the true closest point on the unit sphere (normalized
              // original point)
              const Point<dim> &true_closest_point =
                original_point / original_point.norm();

              // Check that the closest point is actually on the unit sphere
              const double distance_from_origin = closest_point.norm();
              if (std::abs(distance_from_origin - 1.0) < tolerance)
                {
                  n_points_on_sphere++;
                }
              else
                {
                  deallog << "Point not on unit sphere: "
                          << "Original point: " << original_point
                          << ", Closest point: " << closest_point
                          << ", Distance from origin: " << distance_from_origin
                          << std::endl;
                }

              // Check that the closest point is close to the true closest
              // point
              const double distance_to_true =
                (closest_point - true_closest_point).norm();
              if (distance_to_true > tolerance)
                {
                  deallog << "Distance to true closest point too large: "
                          << "Original point: " << original_point
                          << ", Closest point: " << closest_point
                          << ", True closest point: " << true_closest_point
                          << ", Distance: " << distance_to_true << std::endl;
                }

              // Check if the real point corresponds with the unit point
              const Point<dim> unit_closest_point =
                mapping.transform_real_to_unit_cell(cell, closest_point);
              if ((unit_closest_point - closest_unit_reference_points[q])
                    .norm() > 1e-8)
                {
                  deallog << "Closest point does not correspond to unit point: "
                          << "Original point: " << original_point
                          << ", Closest point: " << closest_point
                          << ", Unit closest point: " << unit_closest_point
                          << std::endl;
                }
            }
        }
    }

  deallog << "Tested " << n_tested_points << " points" << std::endl;
  deallog << "Points on unit sphere (within tolerance): " << n_points_on_sphere
          << std::endl;
  deallog << "Success rate: "
          << (double)n_points_on_sphere / n_tested_points * 100.0 << "%"
          << std::endl;
}


int
main()
{
  initlog();

  deallog << std::setprecision(8);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
