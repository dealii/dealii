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



// test the results of FECollection::n_blocks()


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // test things with a collection of
  // primitive elements
  {
    hp::FECollection<dim> fe_collection;
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), dim));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), dim));
    AssertThrow(fe_collection.n_components() == dim, ExcInternalError());
  }

  // now the same with one of the elements
  // being non-primitive. the other one can
  // then not simply be a FESystem but must
  // in fact be an FESystem of FESystem to
  // hide multiple components in one block
  if (dim > 1)
    {
      hp::FECollection<dim> fe_collection;
      fe_collection.push_back(
        FESystem<dim>(FESystem<dim>(FE_Q<dim>(2), dim), 1));
      fe_collection.push_back(FE_RaviartThomas<dim>(1));
      AssertThrow(fe_collection.n_blocks() == 1, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
