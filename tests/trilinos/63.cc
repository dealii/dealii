//----------------------------  trilinos_63.cc  ---------------------------
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
//----------------------------  trilinos_63.cc  ---------------------------


#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::SparseMatrix &m)
{
  Assert (m.m() == 100, ExcInternalError());
  Assert (m.n() == 100, ExcInternalError());

  m = 0;
  
  Assert (m.m() == 100, ExcInternalError());
  Assert (m.n() == 100, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("63/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        TrilinosWrappers::SparseMatrix v (100U,100U,5U);
        test (v);
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
