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



// compute the inverse of a small matrix using the SparseILU with
// infinite fill-in

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>
#include "testmatrix.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/vector.h>


int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  for (unsigned int N=1; N<5; ++N)
    {
      deallog << "N=" << N << std::endl;

      SparsityPattern structure(N, N, N);
      for (unsigned int i=0; i<N; ++i)
        for (unsigned int j=0; j<N; ++j)
          structure.add(i,j);
      structure.compress();
      SparseMatrix<double>  A(structure);
      for (unsigned int i=0; i<N; ++i)
        {
          A.set(i,i,2);
          if (i>=1)
            A.set(i,i-1, -1);
          if (i<N-1)
            A.set(i,i+1, -1);
        }

      SparseILU<double> ilu;
      ilu.initialize (A, SparseILU<double>::AdditionalData());

      // now get an explicit
      // representation of the
      // inverse
      FullMatrix<double> inverse (N,N);
      Vector<double> tmp1(N), tmp2(N);
      for (unsigned int i=0; i<N; ++i)
        {
          tmp1 = 0;
          tmp1(i) = 1;
          ilu.vmult (tmp2, tmp1);
          for (unsigned int j=0; j<N; ++j)
            inverse(i,j) = tmp2(j);
        }

      deallog << "Matrix A:" << std::endl;
      A.print_formatted (deallog.get_file_stream(), 3, false);

      deallog << "Matrix A^{-1}:" << std::endl;
      inverse.print_formatted (deallog.get_file_stream(), 3, false);

      deallog << std::endl;
    }
}

