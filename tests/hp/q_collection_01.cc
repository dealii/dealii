//----------------------------  q_collection_01.cc  ---------------------------
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
//----------------------------  q_collection_01.cc  ---------------------------


// test that QCollection objects are copyable without running into
// troubles when the copy is destroyed earlier than the original
// object


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/hp/q_collection.h>

#include <fstream>


template <int dim>
void test ()
{
  hp::QCollection<dim> q_collection;
  q_collection.push_back (QGauss<dim>(2));
  q_collection.push_back (QGauss<dim>(3));

                                   // now create a copy and make sure
                                   // it goes out of scope before the
                                   // original
  {
    hp::QCollection<dim> copy (q_collection);
  }
}



int main ()
{
  std::ofstream logfile("q_collection_01/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
