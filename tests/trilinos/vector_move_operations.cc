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

#define PRINTME(name, var)                        \
  deallog << "Vector: " name << ':' << std::endl; \
  for (unsigned int i = 0; i < var.size(); ++i)   \
    deallog << var[i] << ' ';                     \
  deallog << std::endl;

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  {
    IndexSet local_owned(10);
    local_owned.add_range(0, 10);

    TrilinosWrappers::MPI::Vector temp(local_owned, MPI_COMM_WORLD);
    for (unsigned int i = 0; i < temp.size(); ++i)
      temp[i] = (double)(i + 1);

    PRINTME("Vector", temp);

    TrilinosWrappers::MPI::Vector u(std::move(temp));
    PRINTME("move constructor", u);

    TrilinosWrappers::MPI::Vector v;
    v = u;
    PRINTME("copy assignment", v);
    PRINTME("old object", u);

    v.clear();
    v = std::move(u);
    PRINTME("move assignment", v);
    deallog << "old object size: " << u.size() << std::endl;
  }
}
