//----------------------------  mapping_collection_04.cc  ---------------------------
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
//----------------------------  mapping_collection_04.cc  ---------------------------


// a test that shows that mapping_collection_0[1-3] really is due to
// the fact that MappingQ has a dysfunctional copy constructor...


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/fe/mapping_q.h>

#include <fstream>


template <int dim>
void test ()
{
  MappingQ<dim> mapping(2);
  {
    deallog << "Copying..." << std::endl;
    MappingQ<dim> copy(mapping);
    deallog << "Deleting clone..." << std::endl;
  }
  deallog << "Destroying original..." << std::endl;  
}



int main ()
{
  std::ofstream logfile("mapping_collection_04/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
