//----------------------------  petsc_64.cc  ---------------------------
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
//----------------------------  petsc_64.cc  ---------------------------


// This test should be run on multiple processors


#include "../tests.h"
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


template<typename MatrixType>
void test (MatrixType &m)
{
  m.add(0,0,1);  
  m = 0;
  m.compress();

  Assert(fabs(m.frobenius_norm())<1e-15, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_64.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
	const unsigned int n_dofs=420;
					 // check
					 // PETScWrappers::SparseMatrix
        PETScWrappers::SparseMatrix
	  v1 (n_dofs, n_dofs, 5);
        test (v1);

					 // check
					 // PETScWrappers::MPI::SparseMatrix	
	MPI_Comm mpi_communicator (MPI_COMM_WORLD);	
	int n_jobs=1;
	MPI_Comm_size (mpi_communicator, &n_jobs);
	const unsigned int n_mpi_processes=static_cast<unsigned int>(n_jobs);
	Assert(n_dofs%n_mpi_processes==0, ExcInternalError());
	const unsigned int n_local_dofs=n_dofs/n_mpi_processes;
        PETScWrappers::MPI::SparseMatrix
	  v2 (mpi_communicator, n_dofs, n_dofs, n_local_dofs, n_local_dofs, 5);
        test (v2);
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
