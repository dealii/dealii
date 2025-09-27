// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that there is at most a one level difference between cells
// across periodic boundary faces for serial Triangulations.

#include <deal.II/base/geometry_info.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  double left  = -1.;
  double right = 1.;

  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria, left, right, true);

  // setup periodic boundary conditions
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;

  // for the hypercube faces
  unsigned int id1 = 0; // id of the left boundary
  unsigned int id2 = 1; // id of the right boundary boundary
  GridTools::collect_periodic_faces(
    tria, id1, id2, /*x-direction*/ 0, periodicity_vector);

  tria.add_periodicity(periodicity_vector);

  const unsigned int n_refinements = 5;
  const double       tolerance     = 1.e-3;

  for (unsigned int n = 0; n < n_refinements; ++n)
    {
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->center()[0] > -tolerance)
          cell->set_refine_flag();

      tria.execute_coarsening_and_refinement();
    }

  // Let's check how many cells have a periodic neighbor.
  unsigned int neighbor_count = 0;

  for (const auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->has_periodic_neighbor(f))
        ++neighbor_count;

  deallog << "Found " << neighbor_count << " faces with periodic neighbor"
          << std::endl;

  const unsigned int expected_neighbor_count =
    Utilities::pow(2, (dim - 1) * (n_refinements - 1)) +
    Utilities::pow(2, (dim - 1) * n_refinements);
  AssertThrow(neighbor_count == expected_neighbor_count, ExcInternalError());
}


int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();
  return 0;
}
