//----------------------------  trilinos_36.cc  ---------------------------
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
//----------------------------  trilinos_36.cc  ---------------------------


// check TrilinosWrappers::Vector::operator-=(Vector) 

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::Vector &v,
           TrilinosWrappers::Vector &w)
{
                                   // set only certain elements of each
                                   // vector
  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      if (i%3 == 0)
        w(i) = i+1.;
    }
  
  v.compress ();
  w.compress ();

  v -= w;

                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      if (i%3 == 0)
        {
          Assert (w(i) == i+1., ExcInternalError());
          Assert (v(i) == i-(i+1.), ExcInternalError());
        }
      else
        {
          Assert (w(i) == 0, ExcInternalError());
          Assert (v(i) == i, ExcInternalError());
        }
    }


  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("36/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);


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
