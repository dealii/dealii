//----------------------------  block_matrix_array_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_matrix_array_02.cc  ---------------------------


// the class BlockMatrixArray had no local type value_type that is
// needed in some places. in particular, this is needed for
// PreconditionBlockSSOR
//
// the test also didn't link before, due to some functions that were
// either missing or in the wrong place
//
// this test doesn't do anything particularly useful, it just makes sure that
// we can compile this kind of code

#include <base/logstream.h>
#include <lac/block_matrix_array.h>
#include <lac/sparse_matrix.h>
#include <lac/precondition_block.h>
#include <iostream>
#include <fstream>


int main () 
{
  std::ofstream logfile("block_matrix_array_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  BlockMatrixArray<SparseMatrix<double> >::value_type i = 1.0;
  deallog << i << std::endl;

				   // the following did not compile
				   // right away
  PreconditionBlockSSOR<BlockMatrixArray<SparseMatrix<double> > > p;
  
  return 0;
}
