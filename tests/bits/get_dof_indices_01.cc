// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// we used to be able to call DoFCellAccessor::get_dof_indices also for
// inactive cells, check that this is now forbidden.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  FE_Q<dim>          fe(1);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  dof_handler.distribute_dofs(fe);

  // loop over all cells, active or
  // not
  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    {
      try
        {
          cell->get_dof_indices(local_dof_indices);
        }
      catch (...)
        {
          deallog << "Assertion: cell not active." << std::endl;
          continue;
        }


      deallog << "Cell = " << cell << ", DoFs=";
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        {
          Assert(local_dof_indices[i] != numbers::invalid_dof_index,
                 ExcInternalError());
          deallog << local_dof_indices[i] << ' ';
        }

      deallog << std::endl;
    }
}



int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
