// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2005 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check that Table<{1,2,3},T>::reinit(TableIndices) work

#include <deal.II/base/table.h>

#include "../tests.h"

template <int dim>
void
test()
{
  dealii::TableIndices<dim> new_size;
  for (unsigned int d = 0; d < dim; ++d)
    new_size[d] = 5;

  dealii::Table<dim, int> unrolled;
  unrolled.reinit(new_size);

  unrolled.clear();
  Assert(unrolled.empty() == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(0);

  test<1>();
  test<2>();
  test<3>();
}
