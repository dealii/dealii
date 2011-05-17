//----------------------------  table_4.cc  ---------------------------
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
//----------------------------  table_4.cc  ---------------------------

// check serialization for Table<4, int>

#include "serialization.h"
#include <deal.II/base/table.h>
#include <boost/serialization/vector.hpp>


void test ()
{ 
  unsigned int index1 = 3, index2 = 4, index3 = 2, index4 = 5;
  TableIndices<4> indices1(index1, index2, index3, index4);
  unsigned int sum_of_indices = index1 + index2 + index3 + index4;
  
  Table<4, int> t1(index1, index2, index3, index4);
  Table<4, int> t2(index1, index2, index3, index4);

  index1 = 2; index2 = 5; index3 = 4; index4 = 1;
  Table<4, int> t3(index1, index2, index3, index4);
  
  unsigned int counter = 0;
  for (unsigned int i1 = 0; i1 < indices1[0]; ++i1)
  {
    for (unsigned int i2 = 0; i2 < indices1[1]; ++i2)
    {
      for (unsigned int i3 = 0; i3 < indices1[2]; ++i3)
      {
        for (unsigned int i4 = 0; i4 < indices1[3]; ++i4)
        {
          t1[i1][i2][i3][i4] = counter ++;
          t2[i1][i2][i3][i4] = counter + sum_of_indices;
        }
      }
    }
  }
  
  verify (t1, t2);
  
  verify (t1, t3);
}


int main ()
{
  std::ofstream logfile("table_4/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
