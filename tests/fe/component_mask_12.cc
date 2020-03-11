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
// here: test conversion from component mask to block mask and back


#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"



void
test()
{
  // the following element has as
  // many components as blocks so
  // that any component mask we come
  // up with should also produce a
  // valid block mask (in fact, it
  // should be exactly the same mask)
  FESystem<2> fe(FE_Q<2>(1), 5);

  // try all possible component
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask = 0; int_mask < (1U << fe.n_components());
       ++int_mask)
    {
      std::vector<bool> component_mask(fe.n_components());
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        component_mask[c] = (int_mask & (1 << c));

      // make sure that the round-trip works
      Assert(ComponentMask(component_mask) ==
               fe.component_mask(fe.block_mask(ComponentMask(component_mask))),
             ExcInternalError());

      // then compare elementwise with
      // the block mask (where there is
      // no direct comparison operator,
      // but they should be the same
      // because each component of the
      // element corresponds to a
      // block)
      for (unsigned int c = 0; c < fe.n_components(); ++c)
        AssertThrow(component_mask[c] == fe.block_mask(component_mask)[c],
                    ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
