// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PreconditionChebyshev on parallel Trilinos matrix and vector


#include <deal.II/base/mpi.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"



int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);


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
      TrilinosWrappers::SparseMatrix AA;
      AA.reinit(A);

      PreconditionChebyshev<TrilinosWrappers::SparseMatrix,
                            TrilinosWrappers::MPI::Vector,
                            TrilinosWrappers::PreconditionJacobi>
        cheby;
      PreconditionChebyshev<
        TrilinosWrappers::SparseMatrix,
        TrilinosWrappers::MPI::Vector,
        TrilinosWrappers::PreconditionJacobi>::AdditionalData cheby_data;
      cheby_data.preconditioner.reset(
        new TrilinosWrappers::PreconditionJacobi());
      cheby_data.preconditioner->initialize(AA);
      cheby_data.degree          = 11;
      cheby_data.smoothing_range = 40;
      cheby.initialize(AA, cheby_data);

      IndexSet set(dim);
      set.add_range(0, dim);
      TrilinosWrappers::MPI::Vector v, tmp1, tmp2;
      v.reinit(set, MPI_COMM_WORLD);
      tmp1.reinit(set, MPI_COMM_WORLD);
      tmp2.reinit(set, MPI_COMM_WORLD);
      for (unsigned int i = 0; i < 3; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            v(j) = random_value<double>();

          AA.vmult(tmp1, v);
          cheby_data.preconditioner->vmult(tmp2, tmp1);
          tmp2 -= v;
          const double ilu_residual = tmp2.l2_norm();

          AA.vmult(tmp1, v);
          cheby.vmult(tmp2, tmp1);
          tmp2 -= v;
          const double cheby_residual = tmp2.l2_norm();

          deallog << "Residual step i=" << i << ":  "
                  << " jacobi=" << ilu_residual << ", cheby=" << cheby_residual
                  << std::endl;
        }
    }

  return 0;
}
