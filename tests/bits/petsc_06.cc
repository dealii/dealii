//----------------------------  petsc_06.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_06.cc  ---------------------------


// check petsc_wrappers::SparseMatrix::l1_norm

#include "../tests.h"
#include <lac/petsc_sparse_matrix.h>    
#include <fstream>
#include <iostream>


void test (petsc_wrappers::SparseMatrix &m)
{
                                   // first set a few entries. count how many
                                   // entries we have
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);  

  m.compress ();

                                   // compare against the exact value of the
                                   // l1-norm (max col-sum)
  deallog << m.l1_norm() << std::endl;
  Assert (m.l1_norm() == 7, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_06.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        petsc_wrappers::SparseMatrix m (5,5,3);
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
