//----------------------------  trilinos_50.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_50.cc  ---------------------------


// check TrilinosWrappers::Vector::operator = (Vector<T>) with T!=TrilinosScalar

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::Vector &v)
{
  Vector<double> w (v.size());
  Vector<float>  x (v.size());

  for (unsigned int i=0; i<w.size(); ++i)
    {
      w(i) = i;
      x(i) = i+1;
    }
  
                                   // first copy from w and make sure we get
                                   // the expected result. then copy from x
                                   // and do the same. in at least one of the
                                   // two cases, the template argument to
                                   // Vector<T> must be different from
                                   // TrilinosScalar
  v = w;
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == i, ExcInternalError());
      Assert (v(i) == i, ExcInternalError());
      Assert (x(i) == i+1, ExcInternalError());
    }

  v = x;
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == i, ExcInternalError());
      Assert (v(i) == i+1, ExcInternalError());
      Assert (x(i) == i+1, ExcInternalError());
    }
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("50/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        TrilinosWrappers::Vector v (100);
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
