//----------------------------  compressed_set_sparsity_pattern_09.cc  ---------------------------
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
//----------------------------  compressed_set_sparsity_pattern_09.cc  ---------------------------


// check CompressedSetSparsityPattern::empty

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <iomanip>
#include <fstream>


void test ()
{
  const unsigned int N = 1000;
  CompressedSetSparsityPattern csp;
  Assert (csp.empty() == true, ExcInternalError());

  csp.reinit (N, N);
  Assert (csp.empty() == false, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("compressed_set_sparsity_pattern_09/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  return 0;
}
