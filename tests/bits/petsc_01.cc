//----------------------------  petsc_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_01.cc  ---------------------------


// check setting elements in a petsc matrix using
// PETScWrappers::SparseMatrix::set()

#include "../tests.h"
#include <lac/petsc_sparse_matrix.h>    
#include <fstream>
#include <iostream>


void test (PETScWrappers::SparseMatrix &m)
{
                                   // first set a few entries
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);

  m.compress ();
  
                                   // then make sure we retrieve the same ones
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          Assert (m(i,j) == i*j*.5+.5, ExcInternalError());
          Assert (m.el(i,j) == i*j*.5+.5, ExcInternalError());
        }
      else
        {
          Assert (m(i,j) == 0, ExcInternalError());
          Assert (m.el(i,j) == 0, ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        PETScWrappers::SparseMatrix m (5,5,3);
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
