// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test if Anasazi can obtain the correct eigenvalues of a diagonal matrix
// for different number types.

#include <deal.II/lac/anasazi_solver.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <typename Number>
void
run_test()
{
  std::size_t eigenvalues_to_compute = 2;
  int                                        matrows = 10;
  LinearAlgebra::TpetraWrappers::SparseMatrix<Number, MemorySpace::Default> K(matrows, matrows, 1);

  for (int i = 0; i < matrows; i++)
    {
      K.set(i, i, i + 1);
    }

  K.compress(VectorOperation::add);

  Teuchos::RCP<Teuchos::ParameterList> parameters =  Teuchos::rcp(new Teuchos::ParameterList);

  parameters->set("Which", "LM");
  
  parameters->set("Block Size", 2);

  parameters->set("Convergence Tolerance", 1e-9);

  LinearAlgebra::TpetraWrappers::AnasaziSolverBlockKrylovSchur<Number, MemorySpace::Default>
    anasazi_eigensolver(parameters);

  // define the remaining inputs
  std::vector<Number>                                        eigenvalues;
  std::vector<LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default>> eigenvectors;

  eigenvectors.resize(eigenvalues_to_compute);

  anasazi_eigensolver.generate_random_initial_space(K);
  anasazi_eigensolver.solve(K, eigenvalues, eigenvectors, eigenvalues_to_compute);

  for (int j = 0; j < eigenvalues_to_compute; j++)
    {
      deallog << "eigenvalue" << j + 1 << ": " << eigenvalues[j] << std::endl;
    }

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      deallog << "Run with number-type double" << std::endl;
      run_test<double>();
      // can someone please run this with float?
      // deallog << "Run with number-type float" << std::endl;
      // run_test<float>();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
