//----------------------------  trilinos_vector_assign_01.cc  ---------------------------
//    $Id: vector_assign_01.cc 24425 2011-09-26 13:53:32Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_vector_assign_01.cc  ---------------------------

/*
  test ::Vector::swap() (fixed in r 25668)
 */


#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>

void print(TrilinosWrappers::Vector &v)
{
  deallog << "size= " << v.size()
	  << " el(0)= " << v(0)
	  << " l2norm()= " << v.l2_norm() << std::endl;
}


void test ()
{
  TrilinosWrappers::Vector v(5);
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = 1;
  TrilinosWrappers::Vector w(9);
  for (unsigned int i=0; i<w.size(); ++i)
    w(i) = 2;
  
  
  deallog << "v: "; print(v);
  deallog << "w: "; print(w);

  deallog << "**swap**" << std::endl;
  
  swap(v,w);

  deallog << "v: "; print(v);
  deallog << "w: "; print(w);

  Assert (v.size()==9, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("vector_swap_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      test ();
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
