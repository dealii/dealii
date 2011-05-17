//----------------------------  fe_collection_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_collection_01.cc  ---------------------------


// test that QCollection objects are copyable without running into
// troubles when the copy is destroyed earlier than the original
// object


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>


template <int dim>
void test ()
{
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (FE_Q<dim>(2));
  fe_collection.push_back (FE_Q<dim>(3));

                                   // now create a copy and make sure
                                   // it goes out of scope before the
                                   // original
  {
    hp::FECollection<dim> copy (fe_collection);
  }
}



int main ()
{
  std::ofstream logfile("fe_collection_01/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
