//----------------------------  block_sparse_matrix_iterator_05.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparse_matrix_iterator_05.cc  ---------------------------


// this tests a failure in the design of the block sparse matrix iterators: falling
// off the end of the matrix does not yield the iterator provided by the end()
// function

#include "../tests.h"
#include <lac/block_sparsity_pattern.h>
#include <lac/block_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test ()
{
                                   // create a 1x2 block matrix with
                                   // non-quadratic blocks so as to make them
                                   // not specially store the diagonal
  BlockSparsityPattern bsp (1,2);
  for (unsigned int j=0; j<2; ++j)
    bsp.block(0,j).reinit (3,2,1);
  bsp.collect_sizes ();

                                   // leave row 0 of block 0,0 empty, but have
                                   // something in this row for block 0,1
  bsp.block(0,0).add (1,0);
  bsp.block(0,1).add (0,0);  
  bsp.compress ();

  BlockSparseMatrix<double> m(bsp);

                                   // get the start iterator. it should point
                                   // to the first element of the first row,
                                   // which happens to be in block 0,1
  BlockSparseMatrix<double>::const_iterator it = m.begin();

  deallog << it->row() << ' '
          << it->column() << ' ' << it->block_row() << ' '
          << it->block_column()
          << std::endl;

  Assert (it->row() == 0, ExcInternalError());
  Assert (it->column() == 2, ExcInternalError());
  Assert (it->block_row() == 0, ExcInternalError());
  Assert (it->block_column() == 1, ExcInternalError());
  
                                   // now advance by two (the only two
                                   // elements of the matrix) and make sure
                                   // that we equal the end iterator
  ++it;
  ++it;
  Assert (it == m.end(), ExcInternalError());
   
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("block_sparse_matrix_iterator_05.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      test ();
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
