//----------------------------  petsc_sparse_matrix_vector_01.cc  ---------------------------
//    $Id: sparse_matrix_vector_01.cc 29164 2013-04-03 16:42:17Z heister $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_sparse_matrix_vector_01.cc  ---------------------------


// document that SparseMatrix::operator= is broken

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  PETScWrappers::SparseMatrix m(10,10,10);
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
        m.set (i,j, i+2*j);


  m.compress();
  
  {
    
  
  PETScWrappers::SparseMatrix m2;

  Mat mm=m;
  Mat m2m=m2;
  
  deallog << &mm << " " << &m2m << std::endl;

  m2 = m;
  {

    Mat mm=m;
    Mat m2m=m2;
  deallog << &mm << " " << &m2m << std::endl;
      
  }
  
  }
  
}


int main (int argc, char **argv)
{
  std::ofstream logfile("sparse_matrix_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test ();
  
	deallog << "OK" << std::endl;
      }
      
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
