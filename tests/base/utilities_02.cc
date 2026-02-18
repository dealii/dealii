// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test Utilities::fixed_power

#include <deal.II/base/utilities.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << Utilities::fixed_power<dim>(2) << std::endl;
  deallog << Utilities::fixed_power<dim>(-2) << std::endl;
  deallog << Utilities::fixed_power<dim>(2.5) << std::endl;
  deallog << Utilities::fixed_power<dim>(-2.5) << std::endl;
  deallog << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
  test<4>();
}
