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


// Check DoFRenumbering::Cuthill_McKee() for FE_DGQ on a refined mesh. The
// sparsity pattern is block diagonal with one connected component per
// cell, which used to result in quadratic complexity.


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(8);

  DoFHandler<2> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE_DGQ<2>(1));
  deallog << "n_dofs: " << dof_handler.n_dofs() << std::endl;

  DoFRenumbering::Cuthill_McKee(dof_handler);

  // every cell has to be numbered consecutively
  std::vector<types::global_dof_index> dof_indices(
    dof_handler.get_fe().n_dofs_per_cell());
  bool consecutive = true;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(dof_indices);
      const auto [min, max] =
        std::minmax_element(dof_indices.begin(), dof_indices.end());
      if (*max - *min + 1 != dof_indices.size())
        consecutive = false;
    }
  deallog << "consecutive: " << std::boolalpha << consecutive << std::endl;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->active_cell_index() == 4)
        break;
      cell->get_dof_indices(dof_indices);
      for (const auto i : dof_indices)
        deallog << i << ' ';
      deallog << std::endl;
    }
}
