//----------------------------  petsc_full_matrix_05.cc  ---------------------------
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
//----------------------------  petsc_full_matrix_05.cc  ---------------------------


// check querying the number of nonzero elements in
// PETScWrappers::FullMatrix

#include "../tests.h"
#include <lac/petsc_full_matrix.h>    
#include <fstream>
#include <iostream>


void test (PETScWrappers::FullMatrix &m)
{
                                   // first set a few entries. count how many
                                   // entries we have
  unsigned int counter = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          m.set (i,j, i*j*.5+.5);
          ++counter;
        }

  m.compress ();

  deallog << m.n_nonzero_elements() << std::endl;
  Assert (m.n_nonzero_elements() == counter,
          ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_full_matrix_05.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        PETScWrappers::FullMatrix m (5,5);
        test (m);
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
