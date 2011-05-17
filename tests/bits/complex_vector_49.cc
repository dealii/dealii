//----------------------------  complex_vector_49.cc  ---------------------------
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
//----------------------------  complex_vector_49.cc  ---------------------------


// check copy constructor Vector<std::complex<double> >::Vector(Vector) for same template
// argument

#include "../tests.h"
#include <deal.II/lac/vector.h>    
#include <fstream>
#include <iomanip>
#include <vector>


void test (Vector<std::complex<double> > &v)
{
                                   // set some entries of the vector
  for (unsigned int i=0; i<v.size(); ++i)
    if (i%3 == 0)
      v(i) = std::complex<double> (i+1., i+2.);
  v.compress ();

                                   // then copy it
  Vector<std::complex<double> > w (v);

                                   // make sure they're equal
  deallog << std::abs(v*w) << ' ' << v.l2_norm() * w.l2_norm()
          << ' ' << std::abs(v*w) - v.l2_norm() * w.l2_norm() << std::endl;
  Assert (std::abs(std::abs(v*w) - v.l2_norm() * w.l2_norm())
	  < 1e-14*std::abs(v*w),
          ExcInternalError());

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("complex_vector_49/output");
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
