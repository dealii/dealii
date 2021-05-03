// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// verify restrictions on level differences imposed by
// hp::Refinement::limit_p_level_difference()
// on h-coarsened grids


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test(const unsigned int max_difference)
{
  Assert(max_difference > 0, ExcInternalError());

  // setup FE collection
  const unsigned int fes_size = 2 * max_difference + 2;

  hp::FECollection<dim> fes;
  while (fes.size() < fes_size)
    fes.push_back(FE_Q<dim>(1));

  const unsigned int contains_fe_index = 0;
  const auto         sequence = fes.get_hierarchy_sequence(contains_fe_index);

  // set up line grid
  // - refine central cell and flag them for coarsening
  // - assign highest p-level to leftmost cell
  //
  // +-------+---+---+-------+
  // |       | c | c |       |
  // |       +---+---+       |
  // |       | c | c |       |
  // +-------+---+---+-------+
  //
  // after prepare_coarsening_and_refinement(), the p-levels should propagate
  // through the central cells as if they were already coarsened

  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 3);

  // refine the central cell
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center()[0] > 1. && cell->center()[0] < 2.)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // now flag these cells for coarsening
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center()[0] > 1. && cell->center()[0] < 2.)
      cell->set_coarsen_flag();

  DoFHandler<dim> dofh(tria);
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->center()[0] < 1.)
      cell->set_active_fe_index(sequence.back());
  dofh.distribute_dofs(fes);

  const bool fe_indices_changed =
    hp::Refinement::limit_p_level_difference(dofh,
                                             max_difference,
                                             contains_fe_index);
  (void)fe_indices_changed;
  Assert(fe_indices_changed, ExcInternalError());

  deallog << "future FE indices before adaptation:" << std::endl;
  for (const auto &cell : dofh.active_cell_iterators())
    deallog << " " << cell->id().to_string() << " " << cell->future_fe_index()
            << std::endl;

  tria.execute_coarsening_and_refinement();

  deallog << "active FE indices after adaptation:" << std::endl;
  for (const auto &cell : dofh.active_cell_iterators())
    deallog << " " << cell->id().to_string() << " " << cell->active_fe_index()
            << std::endl;

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2>(1);
  test<2>(2);
  test<2>(3);
}
