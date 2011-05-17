//----------------------------  compressed_set_sparsity_pattern_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  compressed_set_sparsity_pattern_01.cc  ---------------------------


// check CompressedSetSparsityPattern::n_nonzero_elements. test n_rows() and
// n_cols() in passing

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <iomanip>
#include <fstream>


void test ()
{
                                   // set up a sparsity pattern. since
                                   // CompressedSetSparsityPatterns are most
                                   // often used for 3d, use a rather large
                                   // number of entries per row
  const unsigned int N = 1000;
  CompressedSetSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<40; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);

  Assert (csp.n_rows() == N, ExcInternalError());
  Assert (csp.n_cols() == N, ExcInternalError());
  deallog << csp.n_nonzero_elements () << std::endl;
}



int main () 
{
  std::ofstream logfile("compressed_set_sparsity_pattern_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
