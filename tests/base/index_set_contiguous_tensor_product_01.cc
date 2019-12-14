// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// test IndexSet::contiguous_tensor_product()

#include <deal.II/base/index_set.h>

#include "../tests.h"


int
main()
{
  initlog();
  IndexSet set(10);
  set.add_index(2);
  set.add_range(4, 7);

  set.print(deallog);

  const auto tensor2 = set.tensor_product(complete_index_set(2));
  deallog << "Size: " << tensor2.size() << ": ";
  tensor2.print(deallog);


  const auto tensor3 = set.tensor_product(complete_index_set(3));
  deallog << "Size: " << tensor3.size() << ": ";
  tensor3.print(deallog);

  IndexSet other(5);
  other.add_range(0, 3);
  const auto tensor4 = set.tensor_product(other);
  deallog << "Size: " << tensor4.size() << ": ";
  tensor4.print(deallog);
}
