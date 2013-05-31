//----------------------------  chunk_sparse_matrix_iterator_09.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  chunk_sparse_matrix_iterator_09.cc  ---------------------------


// this test is sparse_matrix_iterator_09 for a ChunkSparseMatrix

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test (const unsigned int chunk_size)
{
  deallog << "Chunk size: " << chunk_size << std::endl;

                                   // create a sparsity pattern with totally
                                   // empty lines (not even diagonals, since
                                   // not quadratic)
  ChunkSparsityPattern sparsity(4,5,1,chunk_size);
  sparsity.add (1,1);
  sparsity.add (3,1);
  sparsity.compress ();

                                   // attach a sparse matrix to it
  ChunkSparseMatrix<double> A(sparsity);

                                   // and loop over the elements of it
  for (ChunkSparseMatrix<double>::const_iterator k=A.begin();
       k!=A.end(); ++k)
    deallog << k->row() << ' ' << k->column() << ' ' << k->value()
            << std::endl;
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_iterator_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test (1);
      test (3);
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
