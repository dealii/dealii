//----------------------------  complex_vector_vector_01.cc  ---------------------------
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
//----------------------------  complex_vector_vector_01.cc  ---------------------------


// check existence Vector<std::complex<double> >::Vector(Vector<std::complex<float> >). this conversion
// constructor was disabled previously altogether because of a compiler defect
// that did not honor the 'explicit' keyword on template constructors. this is
// now autoconf'ed.

#include "../tests.h"
#include <deal.II/lac/vector.h>    
#include <fstream>
#include <iomanip>


void test (Vector<std::complex<double> > &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = std::complex<double> (i+1., i+2.);
  Vector<std::complex<float> > w(v);

  Assert (w==v, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("complex_vector_vector_01/output");
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
