//----------------------------  trilinos_01.cc  ---------------------------
//    $Id: 01.cc 24425 2011-09-26 13:53:32Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_01.cc  ---------------------------


// test contructing vectors

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/abstract_linear_algebra.h>
#include <fstream>
#include <iostream>

template <class VEC>
void test ()
{
  VEC x(3);
  x(0)=10;
  deallog << x.l2_norm() << std::endl;
  x.compress();
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("ala_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      {
        test<dealii::LinearAlgebraPETSc::Vector> ();
        test<dealii::LinearAlgebraDealII::Vector> ();
	test<dealii::LinearAlgebraTrilinos::Vector> ();
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
