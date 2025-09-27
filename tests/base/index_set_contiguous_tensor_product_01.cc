// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
