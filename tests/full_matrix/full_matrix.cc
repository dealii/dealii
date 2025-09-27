// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/eigen.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

const double entries[9] = {11, 12, 13, 21, 22, 23, 31, 32, 33};

// Create a positive definite random matrix

void
random_matrix(FullMatrix<double> &A)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    for (unsigned int j = 0; j < A.n(); ++j)
      {
        double rnd = Testing::rand();
        rnd /= RAND_MAX;
        A(i, j) = (i == j) ? A.m() + rnd : rnd;
      }
}

int
main()
{
  initlog();
  deallog << std::fixed;
  deallog << std::setprecision(3);
  Testing::srand(3391466);

  std::ostream &logfile = deallog.get_file_stream();

  FullMatrix<double> T(3, 3, entries);
  T.print_formatted(logfile, 0, false);

  FullMatrix<double>::const_iterator it  = T.begin();
  FullMatrix<double>::const_iterator end = T.end();
  while (it != end)
    {
      deallog << "Row " << it->row() << "\tCol " << it->column() << "\tVal "
              << it->value() << std::endl;
      ++it;
    }

  it  = T.begin(1);
  end = T.end(1);
  while (it != end)
    {
      deallog << "Row " << it->row() << "\tCol " << it->column() << "\tVal "
              << it->value() << std::endl;
      ++it;
    }

  for (unsigned int i = 1; i < 10; ++i)
    {
      FullMatrix<double> A(i, i), B(i, i);

      // Create matrix and its inverse
      random_matrix(A);
      B.invert(A);

      // Check if unit vectors are recovered
      deallog << "Inverse(dim=" << i << "):";
      for (unsigned int j = 0; j < i; ++j)
        {
          Vector<double> x(i);
          Vector<double> y(i);
          Vector<double> z(i);
          x(j) = 1.;
          A.vmult(y, x);
          B.vmult(z, y);
          z.add(-1., x);
          double a = z.l2_norm();
          if (a > 1.e-12)
            deallog << a << ' ';
        }
      deallog << std::endl;
    }

  if (true)
    {
      FullMatrix<double> A(5, 5), C(5, 5), D(5, 5), H(5, 5);
      D(0, 0) = 1.;
      D(1, 1) = 2.;
      D(2, 2) = 3.;
      D(3, 3) = 4.;
      D(4, 4) = 5.;

      A = D;

      for (unsigned int i = 0; i < 4; ++i)
        {
          // Setup rotation matrix
          C = 0;
          C.diagadd(1.);
          C(i, i) = C(i + 1, i + 1) = std::cos(i + 1.);
          C(i + 1, i)               = std::sin(i + 1.);
          C(i, i + 1)               = -std::sin(i + 1.);

          C.print_formatted(logfile, 3, false);
          deallog << "l1-norm: " << C.l1_norm() << std::endl;
          D = C;
          D.gauss_jordan();
          D.print_formatted(logfile, 3, false);
          deallog << "linfty-norm: " << D.linfty_norm() << std::endl
                  << "Frobenius-norm: " << D.frobenius_norm() << std::endl;

          // Rotate original matrix
          A.mmult(H, C);
          C.Tmmult(A, H);
        }

      A.print_formatted(logfile, 3, false);

      Vector<double>                      u(5);
      GrowingVectorMemory<Vector<double>> mem;

      SolverControl control(500, 1.e-8, false, false);

      if (true)
        {
          u = 1.;
          EigenPower<Vector<double>> von_Mises(control, mem, 0.);
          double                     eigen = 0.;
          von_Mises.solve(eigen, A, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double>> von_Mises(control, mem, -4.);
          double                     eigen = 0.;
          von_Mises.solve(eigen, A, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      H = A;
      H.gauss_jordan();
      H.print_formatted(logfile, 3, false);
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double>> von_Mises(control, mem, 0.);
          double                     eigen = 0.;
          von_Mises.solve(eigen, H, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double>> von_Mises(control, mem, -4.);
          double                     eigen = 0.;
          von_Mises.solve(eigen, H, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
    }
}
