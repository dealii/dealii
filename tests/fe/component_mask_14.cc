// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
      ComponentMask component_mask(fe.n_components(), false);
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        component_mask.set(c, (int_mask & (1 << c)));

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
