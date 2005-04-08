//----------------------------  petsc_sparse_matrix_iterator_01.cc  ---------------------------
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
//----------------------------  petsc_sparse_matrix_iterator_01.cc  ---------------------------


// test PETScWrappers::MatrixBase::const_iterator

#include "../tests.h"
#include <lac/petsc_sparse_matrix.h>    
#include <fstream>
#include <iostream>


void test ()
{
  PETScWrappers::SparseMatrix m(5,5,5);
  m.set (0,0,1);
  m.set (1,1,2);
  m.set (1,2,3);  
  m.compress ();
  PETScWrappers::SparseMatrix::const_iterator i = m.begin();
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  ++i;
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  i++;
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("petsc_sparse_matrix_iterator_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        test ();
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
