//----------------------------  trilinos_49.cc  ---------------------------
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
//----------------------------  trilinos_49.cc  ---------------------------


// like 49, but do the test for
//  TrilinosWrappers::BlockVector
//         ::operator = (dealii::BlockVector<TrilinosScalar>)
// with block vectors instead of plain vectors

#include "../tests.h" 
#include <base/utilities.h>
#include <lac/trilinos_block_vector.h>
#include <lac/block_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::BlockVector &v)
{
  std::vector<unsigned int> sizes (2, 3);
  dealii::BlockVector<TrilinosScalar> w (sizes);

  for (unsigned int i=0; i<w.size(); ++i)
    w(i) = i;
  
  v = w;

  
                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == i, ExcInternalError());
      Assert (v(i) == i, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("49a/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
	std::vector<unsigned int> sizes (2, 3);
        TrilinosWrappers::BlockVector v (sizes);
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
