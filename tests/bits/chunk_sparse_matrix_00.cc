//----------------------------  chunk_sparse_matrix_00.cc  ---------------------------
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
//----------------------------  chunk_sparse_matrix_00.cc  ---------------------------


// set a few elements in a chunk sparse matrix and output them again. should
// yield the same result for all chunk sizes, of course

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>    
#include <fstream>
#include <iomanip>


void test (const unsigned int chunk_size)
{
  deallog << "Chunk size = " << chunk_size << std::endl;
  
  ChunkSparsityPattern sp (5,5,3,chunk_size,false);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  ChunkSparseMatrix<double> m(sp);
  
                                   // first set a few entries
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);

                                   // then write them to the output stream
  for (unsigned int i=0; i<m.m(); ++i)
    {
      for (unsigned int j=0; j<m.n(); ++j)
	deallog << std::setprecision(2) << std::fixed << std::setw(4)
		<< m.el(i,j) << ' ';

      deallog << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_00/output");
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
