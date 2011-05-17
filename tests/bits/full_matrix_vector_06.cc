//----------------------------  full_matrix_vector_06.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_vector_06.cc  ---------------------------


// check FullMatrix::matrix_norm_square

#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<double> &v)
{
  FullMatrix<double> m(v.size(), v.size());
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
        m(i,j) = ( i+2*j);

  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i;
  
  v.compress ();

                                   // <w,Mv>
  const double s = m.matrix_norm_square (v);

                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    Assert (v(i) == i, ExcInternalError());

  double result = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      result += (i+2*j)*j*i;

  Assert (s == result, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("full_matrix_vector_06/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Vector<double> v (100);
      test (v);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
