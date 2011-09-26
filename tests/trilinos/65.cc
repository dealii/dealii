//----------------------------  trilinos_65.cc  ---------------------------
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
//----------------------------  trilinos_65.cc  ---------------------------


// This test used to fail after upgrading to petsc 2.2.1


#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  TrilinosWrappers::MPI::Vector 
    v (Epetra_Map(100,0,Utilities::Trilinos::comm_world()));
  v(0) = 1;
  v = 0;

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("65/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        test ();
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
