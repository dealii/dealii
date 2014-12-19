// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// Test BlockCompressedSimpleSparsityPattern column_number


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main()
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  BlockCompressedSimpleSparsityPattern csp(2,2);
  csp.block(0,0).reinit(2,2);
  csp.block(1,0).reinit(1,2);
  csp.block(0,1).reinit(2,1);
  csp.block(1,1).reinit(1,1);
  csp.collect_sizes();
  for (int i=0; i<3; ++i)
    csp.add(i,i);
  csp.add(1,0);
  csp.add(1,2);
  csp.compress();

  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols() << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1,0).n_rows() << "x" << csp.block(1,0).n_cols() << std::endl;

  csp.print(logfile);

  deallog << std::endl;
  for (unsigned int i=0; i<3; ++i)
    {
      deallog << "row " << i << ": ";

      for (unsigned int j=0; j<csp.row_length(i); ++j)
        deallog << csp.column_number(i,j) << " ";

      deallog << std::endl;
    }

  return 0;
}
