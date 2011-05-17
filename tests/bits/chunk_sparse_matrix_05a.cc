//----------------------------  chunk_sparse_matrix_05a.cc  ---------------------------
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
//----------------------------  chunk_sparse_matrix_05a.cc  ---------------------------


// check querying the number of nonzero elements in
// ChunkSparseMatrix when we don't store the diagonal elements explicitly

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>    
#include <fstream>


void test (const unsigned int chunk_size)
{
  ChunkSparsityPattern sp (5,5,3,chunk_size,false);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  ChunkSparseMatrix<double> m(sp);
  
                                   // first set a few entries. count how many
                                   // entries we have. note that for square
                                   // matrices we also always store the
                                   // diagonal element, except when as above
                                   // we set the special flag for the matrix
                                   // sparsity pattern
  unsigned int counter = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          m.set (i,j, i*j*.5+.5);
          ++counter;
        }

  deallog << m.n_nonzero_elements() << std::endl;

  if (chunk_size == 1)
    Assert (m.n_nonzero_elements() == counter,
	    ExcInternalError())
  else
    Assert (m.n_nonzero_elements() >= counter,
	    ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_05a/output");
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
