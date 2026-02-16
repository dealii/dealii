// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
  deallog << "size: " << vec.locally_owned_size() << std::endl;

  vec = 0.;
  for (unsigned int i = 0; i < vec.locally_owned_size(); ++i)
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
