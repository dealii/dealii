//----------------------------  sparsity_pattern_01_x.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern_01_x.cc  ---------------------------


// check SparsityPattern::copy_from(CompressedSetSparsityPattern)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <iomanip>
#include <fstream>


void test ()
{
  const unsigned int N = 100;
  CompressedSetSparsityPattern csp (N,N);
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
  std::ofstream logfile("sparsity_pattern_01_x/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
