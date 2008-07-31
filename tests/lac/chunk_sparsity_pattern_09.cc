//----------------------------  chunk_sparsity_pattern_09.cc  ---------------------------
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
//----------------------------  chunk_sparsity_pattern_09.cc  ---------------------------


// check ChunkSparsityPattern::copy_from

#include "sparsity_pattern_common.h"

int main ()
{
  std::ofstream logfile("chunk_sparsity_pattern_09/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int chunk_sizes[] = { 1, 2, 4, 5, 7 };
  for (unsigned int i=0; i<sizeof(chunk_sizes)/sizeof(chunk_sizes[0]); ++i)
    {
      chunk_size = chunk_sizes[i];
      copy_from_4<ChunkSparsityPattern> ();
    }
}

  
  
