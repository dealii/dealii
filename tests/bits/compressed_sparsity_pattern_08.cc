//----------------------------  compressed_sparsity_pattern_08.cc  ---------------------------
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
//----------------------------  compressed_sparsity_pattern_08.cc  ---------------------------


// check CompressedSparsityPattern::max_entries_per_row

#include "../tests.h"
#include <base/logstream.h>
#include <lac/compressed_sparsity_pattern.h>
#include <iostream>
#include <fstream>


void test ()
{
                                   // set up a sparsity pattern. since
                                   // CompressedSparsityPatterns are most
                                   // often used for 3d, use a rather large
                                   // number of entries per row
  const unsigned int N = 1000;
  CompressedSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<40; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);

  deallog << ' ' << csp.max_entries_per_row () << std::endl;
}



int main () 
{
  std::ofstream logfile("compressed_sparsity_pattern_08.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
