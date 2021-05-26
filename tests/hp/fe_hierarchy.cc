// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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
