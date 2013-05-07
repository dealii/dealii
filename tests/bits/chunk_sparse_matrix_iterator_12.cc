//----------------------------  chunk_sparse_matrix_iterator_12.cc  ---------------------------
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
//----------------------------  chunk_sparse_matrix_iterator_12.cc  ---------------------------


// like chunk_sparse_matrix_iterator_11 but for const_iterator

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

  ChunkSparseMatrix<double>::const_iterator k = A.begin(),
                                            j = ++A.begin();

  Assert (k < j, ExcInternalError());
  Assert (j > k, ExcInternalError());

  Assert (!(j < k), ExcInternalError());
  Assert (!(k > j), ExcInternalError());

  Assert (k != j, ExcInternalError());
  Assert (!(k == j), ExcInternalError());

  Assert (k == k, ExcInternalError());
  Assert (!(k != k), ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_iterator_12/output");
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
