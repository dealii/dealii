//----------------------------  mapping_collection_02.cc  ---------------------------
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
//----------------------------  mapping_collection_02.cc  ---------------------------


// a test that triggers really hard to track down failures in
// mapping_collection_01 in a really simple way


#include <base/logstream.h>
#include <fe/mapping_collection.h>
#include <fe/mapping_q.h>

#include <fstream>


template <int dim>
void test ()
{
  hp::MappingCollection<dim> mapping_collection;
  mapping_collection.push_back (MappingQ<dim>(2));
}



int main ()
{
  std::ofstream logfile("mapping_collection_02/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
