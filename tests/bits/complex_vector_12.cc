//----------------------------  complex_vector_12.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  complex_vector_12.cc  ---------------------------


// check Vector<std::complex<double> >::operator() in set-mode

#include "../tests.h"
#include <deal.II/lac/vector.h>    
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<std::complex<double> > &v)
{
                                   // set only certain elements of the
                                   // vector. have a bit pattern of where we
                                   // actually wrote elements to
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) = std::complex<double> (i+1., i+2.);
      pattern[i] = true;
    }

  v.compress ();

                                   // check that they are ok, and this time
                                   // all of them
  for (unsigned int i=0; i<v.size(); ++i)
    Assert ((((pattern[i] == true) && (v(i) == std::complex<double> (i+1., i+2.)))
             ||
             ((pattern[i] == false) && (v(i) == std::complex<double> (0)))),
             ExcInternalError());

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("complex_vector_12/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Vector<std::complex<double> > v (100);
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
