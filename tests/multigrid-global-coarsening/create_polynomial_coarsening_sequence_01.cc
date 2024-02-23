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

// Print all polynomial coarsening sequences up to degree 15.

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"


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
