//----------------------------  full_matrix_iterator_01.cc  ---------------------------
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
//----------------------------  full_matrix_iterator_01.cc  ---------------------------


// like sparse_matrix_iterator_12, but for FullMatrix

#include "../tests.h"
#include <lac/full_matrix.h>
#include <fstream>
#include <iostream>


void test ()
{
  FullMatrix<double> A(3,3);

  const FullMatrix<double>::const_iterator k = A.begin(),
                                           j = ++A.begin();

  Assert (k < j, ExcInternalError());
  Assert (j > k, ExcInternalError());

  Assert (!(j < k), ExcInternalError());
  Assert (!(k > j), ExcInternalError());

  Assert (k != j, ExcInternalError());
  Assert (!(k == j), ExcInternalError());

  Assert (k == k, ExcInternalError());
  Assert (!(k != k), ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("full_matrix_iterator_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
