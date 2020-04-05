// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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



// tests for the ComponentMask class
//
// here: test conversion from component mask to block mask does indeed
// not work if we try to split a block


#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"



void
test()
{
  // the following element has more
  // components than blocks
  FESystem<2> fe(FE_RaviartThomas<2>(1), 2);

  // try all possible component
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask = 0; int_mask < (1U << fe.n_components());
       ++int_mask)
    {
      std::vector<bool> component_mask(fe.n_components());
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        component_mask[c] = (int_mask & (1 << c));

      // now try to convert to a block
      // mask. this should not always
      // work
      try
        {
          deallog << fe.block_mask(component_mask) << std::endl;
        }
      catch (ExceptionBase &e)
        {
          deallog << e.get_exc_name() << std::endl;
        }
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(4);

  test();
}
