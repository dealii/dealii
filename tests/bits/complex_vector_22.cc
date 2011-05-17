//----------------------------  complex_vector_22.cc  ---------------------------
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
//----------------------------  complex_vector_22.cc  ---------------------------


// check Vector<std::complex<double> >::operator*(Vector) on two vectors that are
// orthogonal

#include "../tests.h"
#include <deal.II/lac/vector.h>    
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<std::complex<double> > &v,
           Vector<std::complex<double> > &w)
{
                                   // set only certain elements of each
                                   // vector, but disjoint sets of elements
  for (unsigned int i=0; i<v.size(); ++i)
    if (i%3 == 0)
      v(i) = std::complex<double> (i+1., i+2.);
    else
      w(i) = std::complex<double> (i+1., i+2.);
  v.compress ();
  w.compress ();

                                   // make sure the scalar product is zero
  Assert (v*w == std::complex<double>(0), ExcInternalError());

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("complex_vector_22/output");
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
