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
testor(IndexSet &a, IndexSet &other, unsigned int offset, bool verbose)
{
  IndexSet merged(a);

  merged.add_indices(other, offset);

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
      Assert(merged.is_element(i) ==
               (a.is_element(i) ||
                (i >= offset && other.is_element(i - offset))),
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
  testor(id, empty, 2, true);

  deallog << "* add self: " << std::endl;
  testor(id, id, 2, true);

  deallog << "* add id2: " << std::endl;
  IndexSet id2(size);
  id2.add_index(0);
  id2.add_index(2);
  id2.add_index(3);
  testor(id, id2, 3, true);
}



int
main()
{
  initlog();

  test();
}
