// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// Tests generalized eigenvalues of FullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>

const double left[] =
{
  4., -1., -1., -1.,
  -1., 4., -1., -1.,
  -1., -1., 4., -1.,
  -1., -1., -1., 4.
};

const double right[] =
{
  4., -1., -1., -1.,
  -1., 5., -1., -1.,
  -1., -1., 6., -1.,
  -1., -1., -1., 7.
};



int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  logfile.precision(1);
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A(4,4,left),
             B(4,4,right);
  LAPACKFullMatrix<double> LA(4,4),
                   LB(4,4);
  for (unsigned int itype=1; itype<=3; ++itype)
    {
      deallog << std::endl
              << "generalized eigenvalue problem of type "
              << itype << std::endl;
      LA = A;
      LB = B;
      std::vector<Vector<double> > eigenvectors(A.m());
      LA.compute_generalized_eigenvalues_symmetric (LB, eigenvectors, itype);

      for (unsigned int i=0; i<A.m(); ++i)
        {
          std::complex<double> lambda = LA.eigenvalue(i);
          deallog << "generalized eigenvalue "
                  << std::scientific << lambda.real() << '\t'
                  << std::scientific << lambda.imag() << std::endl
                  << "generalized eigenvector ";
          for (unsigned int j=0; j<A.m(); ++j)
            {
              deallog << std::scientific
                      << eigenvectors[i](j)/eigenvectors[i](0)
                      << '\t';
            }
          deallog << std::endl;
        }
    }
}
