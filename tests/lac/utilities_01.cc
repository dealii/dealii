// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test estimate of largest eigenvalue of a matrix.
// The matrix is the same as in slepc/solve_04 test which has eingenvalues:
// 3.98974 > 3.95906 > 3.90828 > 3.83792

#include <deal.II/lac/petsc_compatibility.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/utilities.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <typeinfo>

#include "../tests.h"

#include "../slepc/testmatrix.h"
#include "../testmatrix.h"

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(6);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    const unsigned int size = 31;
    unsigned int       dim  = (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl << std::endl;

    // Make matrix
    FD1DLaplaceMatrix           testproblem(size);
    PETScWrappers::SparseMatrix A(dim, dim, 3);
    testproblem.three_point(A);
    A.compress(VectorOperation::insert);

    PETScWrappers::MPI::Vector v0(MPI_COMM_WORLD, dim, dim);
    PETScWrappers::MPI::Vector y(MPI_COMM_WORLD, dim, dim);
    PETScWrappers::MPI::Vector x(MPI_COMM_WORLD, dim, dim);
    for (unsigned int j = 0; j < v0.size(); ++j)
      v0[j] = random_value<double>();

    v0.compress(VectorOperation::insert);
    GrowingVectorMemory<PETScWrappers::MPI::Vector> vector_memory;

    for (unsigned int k = 4; k < 10; ++k)
      {
        const double est = Utilities::LinearAlgebra::lanczos_largest_eigenvalue(
          A, v0, k, vector_memory);
        Assert(est > 3.98974, ExcInternalError());
        deallog << k << std::endl << "Lanczos " << est << std::endl;

        // estimate from CG
        {
          ReductionControl                     control(k,
                                   std::sqrt(
                                     std::numeric_limits<double>::epsilon()),
                                   1e-10,
                                   false,
                                   false);
          std::vector<double>                  estimated_eigenvalues;
          SolverCG<PETScWrappers::MPI::Vector> solver(control);
          solver.connect_eigenvalues_slot(
            [&estimated_eigenvalues](const std::vector<double> &ev) -> void {
              estimated_eigenvalues = ev;
            });
          y = v0;
          PreconditionIdentity preconditioner;
          try
            {
              solver.solve(A, x, y, preconditioner);
            }
          catch (SolverControl::NoConvergence &)
            {}

          deallog << "CG " << estimated_eigenvalues.back() << std::endl;
        }
      }
  }
}
