//----------------------------  trilinos_vector_assign_02.cc  ---------------------------
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
//----------------------------  trilinos_vector_assign_02.cc  ---------------------------


// this is equivalent to the trilinos_vector_assign_02 test, except that we use
// operator+= instead of operator=. Now, this does not present a problem,
// since the compiler does not automatically generate a version of this
// operator, but simply performs the conversion to TrilinosScalar, i.e. the
// argument to the user-defined operator+=. This is not exciting, but since I
// wrote the test to make sure it works this way, let's keep it then...

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::Vector &v,
           TrilinosWrappers::Vector &w)
{
                                   // set the first vector
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i;

                                   // add elements by reference
  for (unsigned int i=0; i<v.size(); ++i)
    w(i) += v(i);

                                   // check that they're equal
  Assert (v==w, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("vector_assign_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        TrilinosWrappers::Vector v (100);
        TrilinosWrappers::Vector w (100);
        test (v,w);
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
