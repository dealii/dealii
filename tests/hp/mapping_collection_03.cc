//----------------------------  mapping_collection_03.cc  ---------------------------
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
//----------------------------  mapping_collection_03.cc  ---------------------------


// a test that triggers really hard to track down failures in
// mapping_collection_01 in a really simple way; an even further
// simplified form of mapping_collection_02


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/fe/mapping_q.h>

#include <fstream>


template <int dim>
void test ()
{
  MappingQ<dim> mapping(2);
  deallog << "Cloning..." << std::endl;
  Mapping<dim> *copy = mapping.clone();
  deallog << "Deleting clone..." << std::endl;
  delete copy;
  deallog << "Destroying original..." << std::endl;  
}



int main ()
{
  std::ofstream logfile("mapping_collection_03/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
