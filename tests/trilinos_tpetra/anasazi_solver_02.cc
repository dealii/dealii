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

// Test if Anasazi can correctly compute the eigenvalues of a FiniteDifference
// five-point matrix. Basically the same as slepc/solve_01.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/anasazi_solver.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "testmatrix_fd_1d.h"
#include "../testmatrix.h"


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(7);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int size = 46;
  unsigned int       dim  = (size - 1) * (size - 1);

  const int n_eigenvalues = 4;

  deallog << "Size " << size << " Unknowns " << dim << std::endl << std::endl;

  FDMatrix     testproblem(size, size);
  FDDiagMatrix diagonal(size, size);
  LinearAlgebra::TpetraWrappers::SparseMatrix<double> A(
    dim, dim, 5);
  LinearAlgebra::TpetraWrappers::SparseMatrix<double> B(
    dim, dim, 5);

  testproblem.five_point(A);
  A.compress(VectorOperation::insert);
  diagonal.diag(B);
  B.compress(VectorOperation::insert);

  std::vector<
    LinearAlgebra::TpetraWrappers::Vector<double>>
                      u(n_eigenvalues,
      LinearAlgebra::TpetraWrappers::Vector<double>(
        complete_index_set(dim), MPI_COMM_WORLD));
  std::vector<double> v(n_eigenvalues);

  Teuchos::RCP<Teuchos::ParameterList> parameters =  Teuchos::rcp(new Teuchos::ParameterList);

  
  parameters->set("Maximum Restarts", 2000);
  parameters->set("Convergence Tolerance", 1e-6);

  /*int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
  verbosity += Anasazi::IterationDetails;
  parameters->set( "Verbosity", verbosity );*/
  parameters->set("Relative Convergence Tolerance", true);

  // AnasaziSolverBlockKrylovSchur
  {
    parameters->set("Block Size", 1);
    parameters->set("Which", "LM");
    LinearAlgebra::TpetraWrappers::AnasaziSolverBlockKrylovSchur<double> solverkrylov(parameters);
    deallog << "AnasaziSolverBlockKrylovSchur" << std::endl;
    solverkrylov.generate_random_initial_space(A);
    solverkrylov.solve(A, B, v, u, v.size());

    for (int i = 0; i < v.size(); i++)
      {
        deallog << "Eigenvalues" << i + 1 << ": " << v[i] << std::endl;
      }
    
    deallog << "OK" << std::endl;
  }

  // AnasaziSolverGeneralizedDavidson
  {
    parameters->set("Block Size", n_eigenvalues);
    parameters->set("Which", "LM");
    LinearAlgebra::TpetraWrappers::AnasaziSolverGeneralizedDavidson<double>
    solvergeneralizeddavidson(parameters);
    deallog << "AnasaziSolverGeneralizedDavidson" << std::endl;
    solvergeneralizeddavidson.generate_random_initial_space(A);
    solvergeneralizeddavidson.solve(A, B, v, u, v.size());

    for (int i = 0; i < v.size(); i++)
      {
        deallog << "Eigenvalues" << i + 1 << ": " << v[i] << std::endl;
      }
    
    deallog << "OK" << std::endl;
  }

}
