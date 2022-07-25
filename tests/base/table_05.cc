// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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
