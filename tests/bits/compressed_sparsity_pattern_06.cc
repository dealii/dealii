//----------------------------  compressed_sparsity_pattern_06.cc  ---------------------------
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
//----------------------------  compressed_sparsity_pattern_06.cc  ---------------------------


// check CompressedSparsityPattern::print_gnuplot. since we create quite some
// output here, choose smaller number of rows and entries than in the other
// tests

#include "../tests.h"
#include <base/logstream.h>
#include <lac/compressed_sparsity_pattern.h>
#include <iostream>
#include <fstream>


std::ofstream logfile("compressed_sparsity_pattern_06.output");

void test ()
{
  const unsigned int N = 100;
  CompressedSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<10; ++j)
      csp.add (i, (i+(i+1)*(j*j+i))%N);
  csp.symmetrize ();

  csp.print_gnuplot (logfile);
  deallog << "OK" << std::endl;
}



int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ();
  return 0;
}
