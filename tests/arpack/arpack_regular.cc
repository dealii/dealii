// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Test ARPACK regular mode for standard eigenvalue problem

#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

void
test()
{
  const unsigned int size          = 100;
  const unsigned int n_eigenvalues = 10;
  // Build a diagonal matrix and check the eigenvalues and the eigenvectors.
  std::vector<std::vector<unsigned int>> column_indices(
    size, std::vector<unsigned int>(1));
  for (unsigned int i = 0; i < size; ++i)
    column_indices[i][0] = i;
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(size,
                             size,
                             column_indices.begin(),
                             column_indices.end());
  SparseMatrix<double> sparse_matrix(sparsity_pattern);
  for (unsigned int i = 0; i < size; ++i)
    sparse_matrix.diag_element(i) = static_cast<double>(i + 1);

  // Compute the eigenvalues and the eigenvectors
  std::vector<std::complex<double>> eigenvalues(n_eigenvalues);
  std::vector<Vector<double>> eigenvectors(n_eigenvalues, Vector<double>(size));
  SolverControl               solver_control(size, 1e-15);
  const unsigned int          num_arnoldi_vectors = 2 * n_eigenvalues + 2;
  ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors,
                                               ArpackSolver::smallest_magnitude,
                                               true,
                                               ArpackSolver::regular,
                                               ArpackSolver::standard);
  ArpackSolver                 eigensolver(solver_control, additional_data);

  SparseMatrix<double> dummy_1;
  SparseMatrix<double> dummy_2;
  eigensolver.solve(sparse_matrix, dummy_1, dummy_2, eigenvalues, eigenvectors);

  // Compare the result to a reference solution
  std::vector<std::complex<double>> ref_eigenvalues(n_eigenvalues);
  std::vector<Vector<double>>       ref_eigenvectors(n_eigenvalues,
                                                     Vector<double>(size));
  for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      ref_eigenvalues[i]     = static_cast<double>(i + 1);
      ref_eigenvectors[i][i] = 1.;
    }

  const double tol = 1e-13;
  for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      std::cout << eigenvalues[i].real() << " " << ref_eigenvalues[i].real()
                << " " << eigenvalues[i].real() - ref_eigenvalues[i].real()
                << std::endl;
      std::cout << eigenvalues[i].imag() << " " << ref_eigenvalues[i].imag()
                << std::endl;
      AssertThrow(std::abs(eigenvalues[i].real() - ref_eigenvalues[i].real()) <
                    tol,
                  ExcInternalError());
      AssertThrow(std::abs(eigenvalues[i].imag() - ref_eigenvalues[i].imag()) <
                    tol,
                  ExcInternalError());
      for (unsigned int j = 0; j < size; ++j)
        AssertThrow(std::abs(std::abs(eigenvectors[i][j]) -
                             ref_eigenvectors[i][j]) < tol,
                    ExcInternalError());
    }
}

int
main()
{
  initlog();

  test();

  deallog << "OK" << std::endl;
}
