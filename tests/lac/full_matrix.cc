// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/eigen.h>

const double entries[9] = { 11,12,13,21,22,23,31,32,33 };

// Create a positive definite random matrix

void random_matrix(FullMatrix<double> &A)
{
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      {
        double rnd = Testing::rand();
        rnd /= RAND_MAX;
        A(i,j) = (i==j) ? A.m()+rnd : rnd;
      }
}

int
main ()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  Testing::srand(3391466);

  FullMatrix<double> T(3,3,entries);
  T.print_formatted(logfile, 0, false);

  FullMatrix<double>::const_iterator it = T.begin();
  FullMatrix<double>::const_iterator end = T.end();
  while (it != end)
    {
      deallog << "Row " << it->row()
              << "\tCol " << it->column()
              << "\tVal " << it->value()
              << std::endl;
      ++it;
    }

  it = T.begin(1);
  end = T.end(1);
  while (it != end)
    {
      deallog << "Row " << it->row()
              << "\tCol " << it->column()
              << "\tVal " << it->value()
              << std::endl;
      ++it;
    }

  for (unsigned int i=1; i<10; ++i)
    {
      FullMatrix<double> A(i,i), B(i,i);

      // Create matrix and its inverse
      random_matrix(A);
      B.invert(A);

      // Check if unit vectors are recovered
      deallog << "Inverse(dim=" << i <<"):";
      for (unsigned int j=0; j<i; ++j)
        {
          Vector<double> x(i);
          Vector<double> y(i);
          Vector<double> z(i);
          x(j) = 1.;
          A.vmult(y,x);
          B.vmult(z,y);
          z.add(-1.,x);
          double a = z.l2_norm();
          if (a > 1.e-12) deallog << a << ' ';
        }
      deallog << std::endl;
    }

  if (true)
    {
      FullMatrix<double> A(5,5), C(5,5), D(5,5), H(5,5);
      D(0,0) = 1.;
      D(1,1) = 2.;
      D(2,2) = 3.;
      D(3,3) = 4.;
      D(4,4) = 5.;

      A = D;

      for (unsigned int i=0; i<4; ++i)
        {
          // Setup rotation matrix
          C = 0;
          C.diagadd(1.);
          C(i,i) = C(i+1,i+1) = std::cos(i+1.);
          C(i+1,i) = std::sin(i+1.);
          C(i,i+1) = -std::sin(i+1.);

          C.print_formatted (logfile,3,false);
          deallog << "l1-norm: " << C.l1_norm() << std::endl;
          D = C;
          D.gauss_jordan();
          D.print_formatted (logfile,3,false);
          deallog << "linfty-norm: " << D.linfty_norm() << std::endl
                  << "Frobenius-norm: " << D.frobenius_norm() << std::endl;

          // Rotate original matrix
          A.mmult(H,C);
          C.Tmmult(A,H);
        }

      A.print_formatted (logfile,3,false);

      Vector<double> u(5);
      GrowingVectorMemory<Vector<double> > mem;

      SolverControl control (500,1.e-8, false, false);

      if (true)
        {
          u = 1.;
          EigenPower<Vector<double> >
          von_Mises(control, mem, 0.);
          double eigen = 0.;
          von_Mises.solve(eigen, A, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double> >
          von_Mises(control, mem, -4.);
          double eigen = 0.;
          von_Mises.solve(eigen, A, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      H = A;
      H.gauss_jordan();
      H.print_formatted (logfile,3,false);
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double> >
          von_Mises(control, mem, 0.);
          double eigen = 0.;
          von_Mises.solve(eigen, H, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
      if (true)
        {
          u = 1.;
          EigenPower<Vector<double> >
          von_Mises(control, mem, -4.);
          double eigen = 0.;
          von_Mises.solve(eigen, H, u);
          deallog << "Eigenvalue: " << eigen << std::endl;
        }
    }
}


