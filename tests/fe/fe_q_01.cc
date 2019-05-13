// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
