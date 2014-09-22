// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check ChunkSparsityPattern::matrix_position

#include "sparsity_pattern_common.h"

int main ()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int chunk_sizes[] = { 1, 2, 4, 5, 7 };
  for (unsigned int i=0; i<sizeof(chunk_sizes)/sizeof(chunk_sizes[0]); ++i)
    {
      chunk_size = chunk_sizes[i];
      matrix_position<ChunkSparsityPattern> ();
    }
}



