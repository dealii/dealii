// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/dofs/dof_handler.h>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
log_dof_diagnostics(const DoFHandler<dim, spacedim> &dof_handler)
{
  // locally relevant cells and dof indices on each
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (!cell->is_artificial())
      {
        deallog << "Cell: " << cell << std::endl;
        deallog << "  ";
        if (cell->is_locally_owned())
          deallog << "locally owned";
        else if (cell->is_ghost())
          deallog << "ghost";
        else
          Assert(false, ExcInternalError());
        deallog << std::endl;

        std::vector<types::global_dof_index> dof_indices(
          cell->get_fe().dofs_per_cell);
        cell->get_dof_indices(dof_indices);
        deallog << " ";
        for (auto i : dof_indices)
          deallog << ' ' << i;
        deallog << std::endl;
      }

  // dof indices that are locally owned
  deallog << "locally_owned_dofs: ";
  dof_handler.locally_owned_dofs().print(deallog);

  // number of dof indices that are locally owned
  deallog << "n_locally_owned_dofs: " << dof_handler.n_locally_owned_dofs()
          << std::endl;

  // number of all dof indices
  deallog << "n_global_dofs: " << dof_handler.n_dofs() << std::endl;
}
