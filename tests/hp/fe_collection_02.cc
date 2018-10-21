// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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



// test the results of FECollection::n_components()


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
    DEAL_II_AssertThrow(fe_collection.n_components() == dim,
                        ExcInternalError());
  }

  // now the same with one of the elements
  // being non-primitive
  if (dim > 1)
    {
      hp::FECollection<dim> fe_collection;
      fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), dim));
      fe_collection.push_back(FE_RaviartThomas<dim>(1));
      DEAL_II_AssertThrow(fe_collection.n_components() == dim,
                          ExcInternalError());
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
