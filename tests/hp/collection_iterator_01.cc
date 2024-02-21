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



// Loop over entries of FECollection.

#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim = 2;

  hp::FECollection<dim> collection(FE_Q<dim>(1), FE_Q<dim>(2), FE_Q<dim>(3));

  for (const auto &fe : collection)
    deallog << fe.n_dofs_per_cell() << std::endl;

  return 0;
}
