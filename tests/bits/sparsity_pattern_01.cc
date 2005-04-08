//----------------------------  sparsity_pattern_01.cc  ---------------------------
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
//----------------------------  sparsity_pattern_01.cc  ---------------------------


// check SparsityPattern::copy_from(CompressedSparsityPattern)

#include "../tests.h"
#include <base/logstream.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <iostream>
#include <fstream>


void test ()
{
  const unsigned int N = 100;
  CompressedSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<10; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);

  SparsityPattern sp;
  sp.copy_from (csp);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<sp.row_length (i); ++j)
      deallog << i << ' ' << j << ' ' << sp.column_number (i,j)
              << std::endl;  
}



int main () 
{
  std::ofstream logfile("sparsity_pattern_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
