//----------------------------  chunk_sparse_matrix_06.cc  ---------------------------
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
//----------------------------  chunk_sparse_matrix_06.cc  ---------------------------


// check ChunkSparseMatrix::l1_norm

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>    
#include <fstream>


void test (const unsigned int chunk_size)
{
  ChunkSparsityPattern sp (5,5,3,chunk_size);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  ChunkSparseMatrix<double> m(sp);
  
                                   // first set a few entries. count how many
                                   // entries we have
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);  

                                   // compare against the exact value of the
                                   // l1-norm (max col-sum)
  deallog << m.l1_norm() << std::endl;
  Assert (m.l1_norm() == 7, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_06/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      const unsigned int chunk_sizes[] = { 1, 2, 4, 5, 7 };
      for (unsigned int i=0;
	   i<sizeof(chunk_sizes)/sizeof(chunk_sizes[0]);
	   ++i)
	test (chunk_sizes[i]);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
