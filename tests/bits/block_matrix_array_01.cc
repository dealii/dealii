//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------


// This tests the construction of a BlockMatrixArray and outputs the
// entered blocks using print_latex.

#include <base/logstream.h>
#include <lac/block_matrix_array.h>
#include <lac/full_matrix.h>
#include <iostream>
#include <fstream>


int main () 
{
  std::ofstream logfile("block_matrix_array_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A1(4,4);
  FullMatrix<double> A2(4,4);
  FullMatrix<double> B(4,3);
  FullMatrix<double> C(3,3);
  PrimitiveVectorMemory<Vector<double> > mem;
  
  BlockMatrixArray<FullMatrix<double>, double> block(2,2,mem);
  
  block.enter(A1,0,0);
  block.enter(A2,0,0,2,true);
  block.enter(B,0,1,-3.);
  block.enter(B,0,1,-3.,true);
  block.enter(C,1,1,1.,true);

  block.print_latex(deallog);
  
  return 0;
}
