// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test default hierarchy of any hp::FECollection object.
// In the default implementation, the new indices are either the
// succeeding or preceding indices until limits are reached.


#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;

  while (fe_collection.size() < 3)
    {
      // add dummy FE to collection
      fe_collection.push_back(FE_Q<dim>(1));
      deallog << "size:" << fe_collection.size() << std::endl;

      for (unsigned int fe_index = 0; fe_index < fe_collection.size();
           ++fe_index)
        deallog << " idx:" << fe_index
                << " next:" << fe_collection.next_in_hierarchy(fe_index)
                << " prev:" << fe_collection.previous_in_hierarchy(fe_index)
                << std::endl;
    }
}



int
main()
{
  initlog();

  test<1>();

  deallog << "OK" << std::endl;
}
