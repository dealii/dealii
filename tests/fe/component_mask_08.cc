// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
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
// here: ComponentMask::represents_the_all_selected_mask


#include <deal.II/fe/component_mask.h>

#include "../tests.h"



void
test()
{
  // test for an initialized mask
  AssertThrow(ComponentMask(12, false).represents_the_all_selected_mask() ==
                false,
              ExcInternalError());
  // note the semantics of the following as
  // described in the documentation
  AssertThrow(ComponentMask(12, true).represents_the_all_selected_mask() ==
                false,
              ExcInternalError());
  // test for an empty mask
  AssertThrow(ComponentMask().represents_the_all_selected_mask() == true,
              ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
