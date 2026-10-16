/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
 * Copyright (C) 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * ------------------------------------------------------------------------
 */

// Check that periodic refinement and coarsening decisions are matched
// in serial triangulations in two and three dimensions when the flag
// match_refinement_across_periodic_boundaries is set.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
check_periodic_levels(const Triangulation<dim> &triangulation)
{
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face : cell->face_indices())
      if (cell->has_periodic_neighbor(face))
        AssertThrow(cell->level() == cell->periodic_neighbor(face)->level(),
                    ExcMessage("Periodic neighbors have different levels."));
}


template <int dim>
void
test()
{
  using Smoothing = typename Triangulation<dim>::MeshSmoothing;
  Triangulation<dim> triangulation(
    Smoothing(Triangulation<dim>::match_refinement_across_periodic_boundaries));

  GridGenerator::hyper_cube(triangulation, 0., 1., true);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(
    triangulation, 0, 1, /*direction=*/0, periodic_faces);
  triangulation.add_periodicity(periodic_faces);

  triangulation.refine_global(3);

  // Refine only cells next to the left periodic boundary. The matching
  // smoothing flag must refine the corresponding cells on the right.
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face : cell->face_indices())
      if (cell->at_boundary(face) && cell->face(face)->boundary_id() == 0)
        cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  check_periodic_levels(triangulation);

  // Request coarsening only on the right periodic boundary. The decision
  // must be canceled because the cells on the left side are not marked.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->center()[0] > 0.875)
      cell->set_coarsen_flag();

  const types::global_cell_index cells_before_periodic_coarsening =
    triangulation.n_active_cells();
  triangulation.execute_coarsening_and_refinement();

  AssertThrow(triangulation.n_active_cells() ==
                cells_before_periodic_coarsening,
              ExcMessage("Periodic coarsening was not matched."));
  check_periodic_levels(triangulation);

  // Coarsen an interior strip to check that the periodic rule still permits
  // unrelated coarsening.
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->center()[0] > 0.5 && cell->center()[0] < 0.75)
      cell->set_coarsen_flag();

  const types::global_cell_index cells_before_interior_coarsening =
    triangulation.n_active_cells();
  triangulation.execute_coarsening_and_refinement();

  AssertThrow(triangulation.n_active_cells() < cells_before_interior_coarsening,
              ExcMessage("Expected interior cells to coarsen."));
  check_periodic_levels(triangulation);
}


int
main()
{
  initlog();

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
