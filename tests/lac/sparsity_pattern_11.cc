//----------------------------  sparsity_pattern_11.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_pattern_11.cc  ---------------------------


// check SparsityPattern::block_read/write

#include "sparsity_pattern_common.h"

int main ()
{
  std::ofstream logfile("sparsity_pattern_11/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  block_read_write<SparsityPattern> ();
}

  
  
