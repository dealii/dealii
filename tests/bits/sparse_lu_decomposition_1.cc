//----------------------------  sparse_lu_decomposition_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_lu_decomposition_1.cc  ---------------------------

// this file didn't compile at one point in time due to the private
// inheritance of SparseMatrix by SparseLUDecomposition, and the
// associated lack of accessibility of the Subscriptor functions to
// the SmartPointer
//
// it was fixed around 203-05-22


#include "../tests.h"
#include <base/logstream.h>
#include <base/smartpointer.h>
#include <lac/sparse_ilu.h>
#include <fstream>

  

int main () 
{
  std::ofstream logfile("sparse_lu_decomposition_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  SmartPointer<SparseLUDecomposition<double> > sparse_decomp;

  deallog << "OK" << std::endl;
  
  return 0;
}

