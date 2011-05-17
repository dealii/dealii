//----------------------------  complex_vector_41.cc  ---------------------------
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
//----------------------------  complex_vector_41.cc  ---------------------------


// check Vector<std::complex<double> >::sadd(s, Vector) 

#include "../tests.h"
#include <deal.II/lac/vector.h>    
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<std::complex<double> > &v,
           Vector<std::complex<double> > &w)
{
  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      w(i) = std::complex<double> (i+1., i+2.);
    }
  
  v.compress ();
  w.compress ();

  v.sadd (2, w);

                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == std::complex<double> (i+1., i+2.),
	      ExcInternalError());
      Assert (v(i) == 2.*i+std::complex<double> (i+1., i+2.),
	      ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("complex_vector_41/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Vector<std::complex<double> > v (100);
      Vector<std::complex<double> > w (100);
      test (v,w);
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
