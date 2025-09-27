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
// here: ComponentMask::first_selected_component


#include <deal.II/fe/component_mask.h>

#include "../tests.h"



void
test()
{
  std::vector<bool> v(12, false);
  v[3] = true;

  ComponentMask m(v);

  // test for an initialized mask
  Assert(m.first_selected_component() == 3, ExcInternalError());
  Assert(ComponentMask(12, true).first_selected_component() == 0,
         ExcInternalError());
  // test for an empty mask
  Assert(ComponentMask().first_selected_component(12) == 0, ExcInternalError());

  deallog << "OK" << std::endl;

  // the following should yield an exception:
  try
    {
      Assert(ComponentMask(12, true).first_selected_component(13) == 0,
             ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  // as should this:
  try
    {
      Assert(ComponentMask(12, false).first_selected_component() == 0,
             ExcInternalError());
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(4);

  test();
}
