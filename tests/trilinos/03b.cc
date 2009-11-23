//----------------------------  trilinos_03b.cc  ---------------------------
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
//----------------------------  trilinos_03b.cc  ---------------------------


// check setting elements in a petsc matrix using set() and add()
// intermixed. this poses PETSc some problems, since one has to flush some
// buffer in between these two types of operations
//
// in contrast to trilinos_03a, we set and add the same elements here twice, then
// overwrite them again to get the original value back

#include "../tests.h" 
#include <base/utilities.h>
#include <lac/trilinos_sparse_matrix.h>    
#include <fstream>
#include <iostream>


void test (TrilinosWrappers::SparseMatrix &m)
{
                                   // first set a few entries
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);
                                   // then add the same elements again
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.add (i,j, i*j*.5+.5);

  m.compress();

                                   // and overwrite everything again
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);
  
  m.compress ();
  
                                   // then make sure we retrieve the same ones
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          Assert (m(i,j) == i*j*.5+.5, ExcInternalError());
          Assert (m.el(i,j) == i*j*.5+.5, ExcInternalError());
        }
      else
        {
          Assert (m.el(i,j) == 0, ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("03b/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      {
        TrilinosWrappers::SparseMatrix m (5U,5U,3U);
        test (m);
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
