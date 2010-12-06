//----------------------------  table_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  table_1.cc  ---------------------------

// check serialization for Table<1. int>

#include "serialization.h"
#include <base/table.h>
#include <boost/serialization/vector.hpp>

void test ()
{ 
  unsigned int index1 = 3;
  Table<1, int> t1(index1);
  
  Table<1, int> t2(index1);

  unsigned int index3 = 2;
  Table<1, int> t3(index3);

  for (unsigned int i = 0; i< index1; i++)
  {
    t1[i] = i + 1;
    t2[i] = i + 1 + index1;
  }
  verify (t1, t2);
  
  verify (t1, t3);
}


int main ()
{
  std::ofstream logfile("table_1/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
