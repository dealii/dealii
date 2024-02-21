// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Print all possible degree pairs during polynomial transfer.

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"


void
test()
{
  for (unsigned int i = 1; i < 25; ++i)
    {
      for (unsigned int j = 1; j < 25; ++j)
        if (MGTwoLevelTransfer<2, LinearAlgebra::distributed::Vector<double>>::
              fast_polynomial_transfer_supported(i, j))
          deallog << 1 << ' ';
        else
          deallog << "  ";
      deallog << std::endl;
    }
}

int
main()
{
  initlog();

  test();
}
