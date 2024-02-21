// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// document a bug that would fill all memory if you try to create a
// FE_Q<dim>(0).

#include <deal.II/fe/fe_q.h>

#include <iostream>

#include "../tests.h"

int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      FE_Q<3> fe(0);
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}
