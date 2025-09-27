// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet::add_indices(IndexSet)

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"

void
testor(IndexSet &a, IndexSet &other, bool verbose = true)
{
  IndexSet merged(a);

  merged.add_indices(other);

  if (verbose)
    {
      deallog << "Original index set: " << std::endl;
      a.print(deallog);
      deallog << "other index set: " << std::endl;
      other.print(deallog);
      deallog << "merged index set: " << std::endl;
      merged.print(deallog);
    }

  for (unsigned int i = 0; i < merged.size(); ++i)
    {
      Assert(merged.is_element(i) == (a.is_element(i) || other.is_element(i)),
             ExcInternalError());
    }
}



void
test()
{
  const int size = 10;

  IndexSet empty(size);
  IndexSet id(size);

  id.add_index(3);
  id.add_index(4);
  id.add_index(7);

  deallog << "* add empty: " << std::endl;
  testor(id, empty);

  deallog << "* add self: " << std::endl;
  testor(id, id);

  deallog << "* add id2: " << std::endl;
  IndexSet id2(size);
  id2.add_index(0);
  id2.add_index(2);
  id2.add_index(3);
  testor(id, id2);

  deallog << "* random tests... " << std::endl;
  for (unsigned int i = 0; i < 10; ++i)
    {
      const int size = 100;
      IndexSet  a(size);
      IndexSet  b(size);
      for (unsigned int i = 0; i < 9 * a.size() / 10; ++i)
        {
          a.add_index(Testing::rand() % a.size());
          b.add_index(Testing::rand() % a.size());
        }
      testor(a, b, false);
    }
}



int
main()
{
  initlog();

  test();
}
