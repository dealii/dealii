//----------------------------  table_2.cc  ---------------------------
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
//----------------------------  table_2.cc  ---------------------------

// check serialization for Table<2, int>

#include "serialization.h"
#include <base/table.h>
#include <boost/serialization/vector.hpp>


void test ()
{ 
  unsigned int index1 = 3, index2 = 4;
  TableIndices<2> indices1(index1, index2);
  unsigned int sum_of_indices = index1 + index2;
  
  Table<2, int> t1(index1, index2);
  Table<2, int> t2(index1, index2);

  index1 = 2; index2 = 5;
  Table<2, int> t3(index1, index2);
  
  unsigned int counter = 0;
  for (unsigned int i1 = 0; i1 < indices1[0]; ++i1)
  {
    for (unsigned int i2 = 0; i2 < indices1[1]; ++i2)
    {
      t1[i1][i2] = counter ++;
      t2[i1][i2] = counter + sum_of_indices;
    }
  }
  verify (t1, t2);
  
  verify (t1, t3);
}


int main ()
{
  std::ofstream logfile("table_2/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
