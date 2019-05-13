// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Check printing

#include <deal.II/base/index_set.h>

#include <deal.II/lac/read_write_vector.h>

#include "../tests.h"


void
test()
{
  IndexSet is(10);
  is.add_range(1, 3);
  is.add_index(7);
  LinearAlgebra::ReadWriteVector<double> vec(is);
  deallog << "size: " << vec.n_elements() << std::endl;

  vec = 0.;
  for (unsigned int i = 0; i < vec.n_elements(); ++i)
    vec.local_element(i) += i;

  vec.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  test();

  return 0;
}
