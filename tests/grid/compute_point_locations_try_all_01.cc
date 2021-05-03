// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Test GridTools::compute_point_locations_try_all
// Test for corner case: vector of points, some of which
// lie outside the mesh

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/fe_field_function.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test_compute_pt_loc(unsigned int n_points)
{
  deallog << "Testing for dim = " << dim << std::endl;
  deallog << "Testing on: " << n_points << " points." << std::endl;

  // Creating a grid in the square [0,1]x[0,1]
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(std::max(6 - dim, 2));

  // Creating the finite elements needed:
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Creating the random points
  std::vector<Point<dim>> points;

  for (size_t i = 0; i < n_points; ++i)
    {
      points.push_back(random_point<dim>());
      // One every three points is
      // shifted outside the mesh
      if (points.size() % 3 == dim - 1)
        {
          points.back()[0] += 5.0;
        }
    }

  // Initializing the cache
  GridTools::Cache<dim, dim> cache(tria);

  auto cell_qpoint_map =
    GridTools::compute_point_locations_try_all(cache, points);
  const auto &cells   = std::get<0>(cell_qpoint_map);
  const auto &qpoints = std::get<1>(cell_qpoint_map);
  const auto &maps    = std::get<2>(cell_qpoint_map);
  const auto &other_p = std::get<3>(cell_qpoint_map);
  size_t      n_cells = cells.size();

  deallog << "Points found in " << n_cells << " cells" << std::endl;

  unsigned int checked_points = 0;
  // testing if the transformation is correct:
  for (unsigned int i = 0; i < n_cells; ++i)
    {
      const auto &cell      = cells[i];
      const auto &quad      = qpoints[i];
      const auto &local_map = maps[i];

      // Given the std::get<1>(cell_qpoint_map) of the current cell, compute the
      // real points
      FEValues<dim> fev(fe, quad, update_quadrature_points);
      fev.reinit(cell);
      const auto &real_quad = fev.get_quadrature_points();

      for (unsigned int q = 0; q < real_quad.size(); ++q)
        {
          // Check if points are the same as real points
          if (real_quad[q].distance(points[local_map[q]]) > 1e-10)
            deallog << "Error on cell : " << cell << " at local point " << i
                    << ", corresponding to real point " << points[local_map[q]]
                    << ", that got transformed to " << real_quad[q]
                    << " instead." << std::endl;
          else
            ++checked_points;
        }
    }
  for (auto i : other_p)
    if (i % 3 != dim - 2)
      {
        deallog << "Error: point index " << i << " was not found!" << std::endl;
      }

  if (checked_points + other_p.size() == points.size())
    deallog << "All points where processed correctly" << std::endl;
  else
    {
      deallog << "Error: some points where not processed correctly"
              << std::endl;
      deallog << "Points found inside the mesh: " << checked_points
              << std::endl;
      deallog << "Discarded points: " << points.size() - checked_points
              << std::endl;
      deallog << "Points outside the mesh: " << other_p.size() << std::endl;
      deallog << "Total points processed: " << points.size() << std::endl;
    }
  deallog << "Test finished" << std::endl;
}

int
main()
{
  initlog();

  deallog << "Deal.II compute_pt_loc_try_all:" << std::endl;
  test_compute_pt_loc<2>(101);
  test_compute_pt_loc<3>(202);
}
