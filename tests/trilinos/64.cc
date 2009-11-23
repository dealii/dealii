//----------------------------  trilinos_64.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_64.cc  ---------------------------


// This test should be run on multiple processors.


#include "../tests.h" 
#include <base/utilities.h>
#include <lac/trilinos_sparse_matrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


template<typename MatrixType>
void test (MatrixType &m)
{
  m.set(0,0,1.);  
  m.compress();
  m = 0;
  m.compress();

  Assert(fabs(m.frobenius_norm())<1e-15, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("64/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
	const unsigned int n_dofs=420;
					 // check
					 // TrilinosWrappers::SparseMatrix
        TrilinosWrappers::SparseMatrix
	  v1 (n_dofs, n_dofs, 5U);
        test (v1);

					 // check
					 // TrilinosWrappers::SparseMatrix
	const unsigned int n_jobs = 
	  Utilities::Trilinos::get_n_mpi_processes(Utilities::Trilinos::comm_world());
	Assert(n_dofs%n_jobs==0, ExcInternalError());
	const unsigned int n_local_dofs=n_dofs/n_jobs;
	Epetra_Map map(n_dofs, n_local_dofs, Utilities::Trilinos::comm_world());
        TrilinosWrappers::SparseMatrix v2 (map, 5);
        test (v2);
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
