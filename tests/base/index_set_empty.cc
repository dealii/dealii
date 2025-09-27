// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the invariant that IndexSet::is_empty() is the same as
// IndexSet::n_elements()==0.

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(20);

  if (index_set.is_empty() == (index_set.n_elements() == 0))
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Not OK" << std::endl;
    }

  index_set.add_range(2, 4);
  if (index_set.is_empty() == (index_set.n_elements() == 0))
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Not OK" << std::endl;
    }

  index_set.add_index(6);
  if (index_set.is_empty() == (index_set.n_elements() == 0))
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Not OK" << std::endl;
    }

  index_set.pop_back();
  index_set.pop_back();
  index_set.pop_back();

  if (index_set.is_empty() == (index_set.n_elements() == 0))
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "Not OK" << std::endl;
    }
}


int
main()
{
  initlog();

  test();
}
