// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test the results of FECollection::n_blocks(). test the case where elements
// have the same number of components but different numbers of blocks. this
// needs to lead to an assertion.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

#include <fstream>


template <int dim>
void test ()
{
  // now the same with one of the elements
  // being non-primitive. the other one can
  // then not simply be a FESystem but must
  // in fact be an FESystem of FESystem to
  // hide multiple components in one
  // block. this test tries to see what
  // happens if one doesn't do this
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (FE_RaviartThomas<dim>(1));
  fe_collection.push_back (FESystem<dim>(FE_Q<dim>(2),dim));

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



int main ()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
