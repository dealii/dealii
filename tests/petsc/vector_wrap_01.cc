//----------------------------  petsc_vector_assign_01.cc  ---------------------------
//    $Id: vector_assign_01.cc 17443 2008-10-31 19:25:47Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_vector_assign_01.cc  ---------------------------


// Test the constructor PETScWrappers::Vector(const Vec &) that takes an
// existing PETSc vector.

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (PETScWrappers::Vector &v,
           PETScWrappers::Vector &w)
{
                                   // set the first vector
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i;

                                   // copy elements by reference
  for (unsigned int i=0; i<v.size(); ++i)
    w(i) = v(i);

                                   // check that they're equal
  Assert (v==w, ExcInternalError());
  
  v=w;
  
  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("vector_wrap_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      Vec vpetsc;
      int ierr = VecCreateSeq (PETSC_COMM_SELF, 100, &vpetsc);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      {	
        PETScWrappers::Vector v (vpetsc);
        PETScWrappers::Vector w (100);
        test (v,w);
      }
      ierr = VecDestroy(vpetsc);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

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
