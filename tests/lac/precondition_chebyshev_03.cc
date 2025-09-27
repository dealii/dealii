// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PreconditionChebyshev on more complex matrix and preconditioner


#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"



int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);


  for (unsigned int size = 4; size <= 16; size *= 2)
    {
      unsigned int dim = (size - 1) * (size - 1);

      deallog << "Size " << size << " Unknowns " << dim << std::endl;

      // Make matrix
      FDMatrix        testproblem(size, size);
      SparsityPattern structure(dim, dim, 5);
      testproblem.five_point_structure(structure);
      structure.compress();
      SparseMatrix<double> A(structure);
      testproblem.five_point(A);

      PreconditionChebyshev<SparseMatrix<double>,
                            Vector<double>,
                            SparseILU<double>>
                                                               cheby;
      PreconditionChebyshev<SparseMatrix<double>,
                            Vector<double>,
                            SparseILU<double>>::AdditionalData cheby_data;
      cheby_data.preconditioner.reset(new SparseILU<double>());
      cheby_data.preconditioner->initialize(A);
      cheby_data.degree          = 11;
      cheby_data.smoothing_range = 40;
      cheby.initialize(A, cheby_data);

      Vector<double> v(dim);
      Vector<double> tmp1(dim), tmp2(dim);
      for (unsigned int i = 0; i < 3; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            v(j) = random_value<double>();

          A.vmult(tmp1, v);
          cheby_data.preconditioner->vmult(tmp2, tmp1);
          tmp2 -= v;
          const double ilu_residual = tmp2.l2_norm();

          A.vmult(tmp1, v);
          cheby.vmult(tmp2, tmp1);
          tmp2 -= v;
          const double cheby_residual = tmp2.l2_norm();

          deallog << "Residual step i=" << i << ":  "
                  << " ilu=" << ilu_residual << ", cheby=" << cheby_residual
                  << std::endl;
        }
    }

  return 0;
}
