// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// hp::Refinement::limit_p_level_difference() used to ignore updating future
// FE indices when they coincide with the active FE index of a cell.
// This posed a problem with p-coarsening.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  // setup FE collection
  hp::FECollection<dim> fes;
  while (fes.size() < 3)
    fes.push_back(FE_Q<dim>(1));

  // setup two cell triangulation
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fes);

  // set both cells to the same active FE index,
  // and assign a future FE index one higher and one lower
  for (const auto &cell : dofh.active_cell_iterators())
    {
      cell->set_active_fe_index(1);

      if (cell->id().to_string() == "0_0:")
        cell->set_future_fe_index(2);
      else if (cell->id().to_string() == "1_0:")
        cell->set_future_fe_index(0);
      else
        Assert(false, ExcInternalError());
    }

  const bool future_fe_indices_changed =
    hp::Refinement::limit_p_level_difference(dofh);

  (void)future_fe_indices_changed;
  Assert(future_fe_indices_changed, ExcInternalError());

  // display FE indices for all cells
  for (const auto &cell : dofh.active_cell_iterators())
    deallog << cell->id().to_string() << ": active:" << cell->active_fe_index()
            << " future:" << cell->future_fe_index()
            << " flag:" << cell->future_fe_index_set() << std::endl;

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
