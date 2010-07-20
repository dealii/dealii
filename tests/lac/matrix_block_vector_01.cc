//----------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include <base/logstream.h>
#include <lac/matrix_block.h>
#include <lac/sparse_matrix.h>
#include <lac/full_matrix.h>
#include <lac/block_sparsity_pattern.h>

#include <iostream>
#include <fstream>

using namespace dealii;

void test_sparse()
{
  MatrixBlockVector<SparseMatrix<float> > v;
//  v.add(0, 0, "A", 
}


int main(int argc, char** argv)
{
  std::cerr << argv[0] << std::endl;
  
  std::ofstream logfile("matrices/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  
}
