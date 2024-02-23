// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the results of FECollection::n_blocks(). test the case where elements
// have the same number of components but different numbers of blocks. this
// needs to lead to an assertion.


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // now the same with one of the elements
  // being non-primitive. the other one can
  // then not simply be a FESystem but must
  // in fact be an FESystem of FESystem to
  // hide multiple components in one
  // block. this test tries to see what
  // happens if one doesn't do this
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_RaviartThomas<dim>(1));
  fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), dim));

  // we will get an assertion failure in
  // n_blocks here.
  try
    {
      fe_collection.n_blocks();
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog.get_file_stream().precision(2);

  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
