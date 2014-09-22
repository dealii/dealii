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


// Tests eigenvalues of FullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <float.h>
#include <fstream>
#include <iostream>
#include <vector>

const double left[] =
{

  1.75, -0.433012701892219, 0.0, 0.0,
  -0.433012701892219, 1.25, 0.0, 0.0,
  0.0, 0.0, 3.5, -0.5,
  0.0, 0.0, -0.5, 3.5
};



int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  logfile.precision(1);
  deallog.attach(logfile);
  deallog.depth_console(0);

  FullMatrix<double> A(4,4,left);
  LAPACKFullMatrix<double> LA(4,4);
  LA = A;
  FullMatrix<double> eigenvectors;
  Vector<double> eigenvalues(0);

  LA.compute_eigenvalues_symmetric (0.5, 2.5,
                                    2.0*DBL_MIN,
                                    eigenvalues,
                                    eigenvectors);

  for (unsigned int i=0; i<eigenvalues.size(); ++i)
    {
      deallog << "eigenvalue "
              << std::scientific << eigenvalues(i) << std::endl
              << "eigenvector ";
      for (unsigned int j=0; j<A.m(); ++j)
        {
          deallog << std::scientific
                  << eigenvectors(j,i)/eigenvectors(0,i)
                  << '\t';
        }
      deallog << std::endl;
    }
}
