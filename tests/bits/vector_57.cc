//----------------------------  vector_57.cc  ---------------------------
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
//----------------------------  vector_57.cc  ---------------------------


// check Vector::is_non_zero

#include "../tests.h"
#include <lac/vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (Vector<double> &v)
{
                                   // set only certain elements of the
                                   // vector. they are all positive
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) += i;
      pattern[i] = true;
    }

  v.compress ();

                                   // check that the vector is really
                                   // non-negative
  Assert (v.is_non_negative() == true, ExcInternalError());

                                   // then set a single element to a negative
                                   // value and check again
  v(v.size()/2) = -1;
  Assert (v.is_non_negative() == false, ExcInternalError());  

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("vector_57.output");
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
