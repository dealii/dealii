// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PreconditionChebyshev with power iteration rather than SolverCG to
// estimate the eigenvalues (less accurate)


#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"



int
main()
{
  initlog();
  deallog << std::setprecision(10);

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

      using Chebyshev = PreconditionChebyshev<SparseMatrix<double>,
                                              Vector<double>,
                                              SparseILU<double>>;
      Chebyshev                 cheby;
      Chebyshev::AdditionalData cheby_data;
      cheby_data.preconditioner.reset(new SparseILU<double>());
      cheby_data.preconditioner->initialize(A);
      cheby_data.degree          = 11;
      cheby_data.smoothing_range = 40;
      cheby_data.eigenvalue_algorithm =
        Chebyshev::AdditionalData::EigenvalueAlgorithm::power_iteration;
      cheby.initialize(A, cheby_data);

      Chebyshev cheby_4;
      cheby_data.polynomial_type =
        Chebyshev::AdditionalData::PolynomialType::fourth_kind;
      cheby_4.initialize(A, cheby_data);

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

          A.vmult(tmp1, v);
          cheby_4.vmult(tmp2, tmp1);
          tmp2 -= v;
          const double cheby_4_residual = tmp2.l2_norm();

          deallog << "Residual step i=" << i << ":  "
                  << " ilu=" << ilu_residual << ", cheby=" << cheby_residual
                  << ", 4th-kind cheby=" << cheby_4_residual << std::endl;
        }
    }

  return 0;
}
