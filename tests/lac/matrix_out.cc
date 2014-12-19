// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/matrix_out.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>
#include <iomanip>

int main ()
{
  std::ofstream logfile("output");
  logfile << std::fixed;
  logfile << std::setprecision(2);
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // test for a square full matrix
  if (true)
    {
      FullMatrix<double> full_matrix(4,4);
      for (unsigned int i=0; i<4; ++i)
        full_matrix(i,i) = 1;

      MatrixOut matrix_out;
      matrix_out.build_patches (full_matrix, "full_matrix");
      matrix_out.write_gnuplot (logfile);
    };

  // test for a rectangular sparse
  // matrix
  if (true)
    {
      SparsityPattern sparsity (4,8,7);
      for (unsigned int i=0; i<4; ++i)
        for (unsigned int j=0; j<8; ++j)
          if (i!=j)
            sparsity.add (i,j);
      sparsity.compress ();

      SparseMatrix<double> sparse_matrix(sparsity);
      for (unsigned int i=0; i<4; ++i)
        for (unsigned int j=0; j<8; ++j)
          sparse_matrix.set(i,j, static_cast<signed int>(i-j));

      MatrixOut matrix_out;
      matrix_out.build_patches (sparse_matrix, "sparse_matrix",
                                MatrixOut::Options (true));
      matrix_out.write_eps (logfile);
    };

  // test collation of elements
  if (true)
    {
      FullMatrix<double> full_matrix(20,20);
      for (unsigned int i=0; i<20; ++i)
        for (unsigned int j=0; j<20; ++j)
          full_matrix(i,j) = (1.*i*i/20/20-1.*j*j*j/20/20/20);

      MatrixOut matrix_out;
      matrix_out.build_patches (full_matrix, "collated_matrix",
                                MatrixOut::Options (false, 4));
      matrix_out.write_gmv (logfile);
    };
}
