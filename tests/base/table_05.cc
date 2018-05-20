// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check that Table<{1,2,3},T>::reinit(TableIndices) work

#include "../tests.h"

#include <deal.II/base/table.h>

template <int dim>
void
test()
{
  dealii::TableIndices<dim> new_size;
  for(unsigned int d = 0; d < dim; d++)
    new_size[d] = 5;

  dealii::Table<dim, int> unrolled;
  unrolled.reinit(new_size);

  deallog << "OK" << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);

  test<1>();
  test<2>();
  test<3>();
}
