// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test DoFTools::make_cell_patches with parallel::distributed::Triangulation


#include "block_list.h"

template <int dim>
void
test_block_list(const parallel::distributed::Triangulation<dim> &tr,
                const FiniteElement<dim>                        &fe)
{
  deallog << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  for (unsigned int level = 0; level < tr.n_global_levels(); ++level)
    {
      SparsityPattern bl;
      DoFTools::make_cell_patches(bl, dof, level);
      bl.compress();

      for (unsigned int i = 0; i < bl.n_rows(); ++i)
        {
          deallog << "Level " << level << " Block " << std::setw(3) << i;
          std::vector<unsigned int> entries;
          for (SparsityPattern::iterator b = bl.begin(i); b != bl.end(i); ++b)
            entries.push_back(b->column());

          std::sort(entries.begin(), entries.end());

          for (unsigned int i = 0; i < entries.size(); ++i)
            deallog << ' ' << std::setw(4) << entries[i];
          deallog << std::endl;
        }
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;
  deallog.push("2D");
  test_global_refinement_parallel<2>(&test_block_list<2>);
  deallog.pop();
  deallog.push("3D");
  test_global_refinement_parallel<3>(&test_block_list<3>);
  deallog.pop();
}
