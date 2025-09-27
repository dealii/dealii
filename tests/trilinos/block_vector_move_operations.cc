// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"

#define PRINTBLOCK(name, var)                                \
  deallog << "Block vector: " name << ':' << std::endl;      \
  for (unsigned int i = 0; i < var.n_blocks(); ++i)          \
    {                                                        \
      deallog << "[block " << i << " ]  ";                   \
      for (unsigned int j = 0; j < var.block(i).size(); ++j) \
        deallog << var.block(i)[j] << ' ';                   \
      deallog << std::endl;                                  \
    }


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  {
    std::vector<IndexSet> local_owned(5);
    for (auto &index : local_owned)
      {
        index.set_size(2);
        index.add_range(0, 2);
      }

    TrilinosWrappers::MPI::BlockVector temp(local_owned, MPI_COMM_WORLD);
    for (unsigned int i = 0; i < 5; ++i)
      for (unsigned int j = 0; j < 2; ++j)
        temp.block(i)[j] = (double)(10 * i + j);

    PRINTBLOCK("BlockVector", temp);

    TrilinosWrappers::MPI::BlockVector u(std::move(temp));
    PRINTBLOCK("move constructor", u);

    TrilinosWrappers::MPI::BlockVector v;
    v = u;
    PRINTBLOCK("copy assignment", v);
    PRINTBLOCK("old object", u);

    v.reinit(0);
    v = std::move(u);
    PRINTBLOCK("move assignment", v);
    deallog << "old object size: " << u.n_blocks() << std::endl;
  }
}
