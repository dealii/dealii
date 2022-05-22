// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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

// Print all possible degree pairs during polynomial transfer.

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

using namespace dealii;

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
