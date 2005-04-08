//----------------------------  vector_49.cc  ---------------------------
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
//----------------------------  vector_49.cc  ---------------------------


// check copy constructor Vector<double>::Vector(Vector) for same template
// argument

#include "../tests.h"
#include <lac/vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (Vector<double> &v)
{
                                   // set some entries of the vector
  for (unsigned int i=0; i<v.size(); ++i)
    if (i%3 == 0)
      v(i) = i+1.;
  v.compress ();

                                   // then copy it
  Vector<double> w (v);

                                   // make sure they're equal
  deallog << v*w << ' ' << v.l2_norm() * w.l2_norm()
          << ' ' << v*w - v.l2_norm() * w.l2_norm() << std::endl;
  Assert (std::fabs(v*w - v.l2_norm() * w.l2_norm()) < 1e-14*(v*w),
          ExcInternalError());

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("vector_49.output");
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
