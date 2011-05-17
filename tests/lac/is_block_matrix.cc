//----------------------------  block_matrices.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_matrices.cc  ---------------------------


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>

#include <fstream>
#include <iomanip>
#include <algorithm>




void test () 
{
  std::ofstream logfile("is_block_matrix/output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << IsBlockMatrix<SparseMatrix<double> >::value << ' '
	  << IsBlockMatrix<SparseMatrix<float> >::value << ' '
	  << IsBlockMatrix<SparseMatrixEZ<double> >::value << ' '
    	  << IsBlockMatrix<SparseMatrixEZ<float> >::value << std::endl;
  
  deallog << IsBlockMatrix<BlockSparseMatrix<double> >::value << ' '
	  << IsBlockMatrix<BlockSparseMatrix<float> >::value << ' '
	  << IsBlockMatrix<BlockSparseMatrixEZ<double> >::value << ' '
    	  << IsBlockMatrix<BlockSparseMatrixEZ<float> >::value << std::endl;

  deallog << IsBlockMatrix<SparsityPattern>::value << ' '
	  << IsBlockMatrix<CompressedSparsityPattern>::value << ' '
	  << IsBlockMatrix<CompressedSetSparsityPattern>::value << ' '
    	  << IsBlockMatrix<CompressedSimpleSparsityPattern>::value << std::endl;

  deallog << IsBlockMatrix<BlockSparsityPattern>::value << ' '
	  << IsBlockMatrix<BlockCompressedSparsityPattern>::value << ' '
	  << IsBlockMatrix<BlockCompressedSetSparsityPattern>::value << ' '
    	  << IsBlockMatrix<BlockCompressedSimpleSparsityPattern>::value << std::endl;
}




int main ()
{
  try
    {
      test ();
    }
  catch (std::exception &e)
    {
      std::cerr << std::endl << std::endl
	   << "----------------------------------------------------"
	   << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
	   << "Aborting!" << std::endl
	   << "----------------------------------------------------"
	   << std::endl;
				       // abort
      return 2;
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
				       // abort
      return 3;
    };
  
  
  return 0;
}
