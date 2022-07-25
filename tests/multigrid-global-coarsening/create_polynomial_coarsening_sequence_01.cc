// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Print all polynomial coarsening sequences up to degree 15.

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

using namespace dealii;

void
test(
  const MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType type)
{
  for (unsigned int i = 1; i <= 15; ++i)
    {
      const auto sequence =
        MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
          i, type);

      for (const auto i : sequence)
        deallog << i << ' ';
      deallog << std::endl;
    }
  deallog << std::endl;
}

int
main()
{
  initlog();

  test(
    MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::bisect);
  test(MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::
         decrease_by_one);
  test(MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::
         go_to_one);
}
