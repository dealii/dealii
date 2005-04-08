//----------------------------  petsc_60.cc  ---------------------------
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
//----------------------------  petsc_60.cc  ---------------------------


// check ::Vector::operator = (const PETScWrappers::Vector &)

#include "../tests.h"
#include <lac/petsc_vector.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


void test (PETScWrappers::Vector &v)
{
                                   // set only certain elements of the
                                   // vector.
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) += i;
      pattern[i] = true;
    }

  v.compress ();

  Vector<double> w; w=v;
  Vector<float>  x; x=v;
  
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (v(i) == w(i), ExcInternalError());
      Assert (v(i) == x(i), ExcInternalError());
    }
      
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_60.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        PETScWrappers::Vector v (100);
        test (v);
      }
      PetscFinalize();
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
