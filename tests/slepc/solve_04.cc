// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the SLEPc solvers with HEP
// the unit tests is a mirror of
// SLEPc-3.6.1/src/eps/examples/tests/test4.c


#include <deal.II/lac/petsc_compatibility.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <typeinfo>

#include "../tests.h"

#include "../testmatrix.h"
#include "testmatrix.h"

template <typename SolverType, typename MatrixType, typename VectorType>
void
check_solve(SolverType               &solver,
            const unsigned int        solver_number,
            const SolverControl      &solver_control,
            const MatrixType         &A,
            std::vector<VectorType>  &u,
            std::vector<PetscScalar> &v)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  try
    {
      solver.set_problem_type(EPS_HEP);
      // reset vectors and set them as initial space
      // to avoid dependency on random numbers:
      for (unsigned int i = 0; i < u.size(); ++i)
        for (unsigned int j = 0; j < u[i].size(); ++j)
          u[i][j] = random_value<double>();

      for (auto &vector : u)
        vector.compress(VectorOperation::insert);

      solver.set_initial_space(u);
      // solve
      solver.solve(A, v, u, v.size());
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }

  switch (solver_number)
    {
      case 1:
        check_solver_within_range((void)true, solver_control.last_step(), 5, 5);
        break;
      case 2:
        check_solver_within_range((void)true,
                                  solver_control.last_step(),
                                  24,
                                  24);
        break;
      case 3:
        check_solver_within_range((void)true,
                                  solver_control.last_step(),
                                  21,
                                  21);
        break;
      // below the spread is bigger since Slepc 3.8:
      case 4:
        check_solver_within_range((void)true,
                                  solver_control.last_step(),
                                  119,
                                  128);
        break;
      case 5:
        check_solver_within_range((void)true,
                                  solver_control.last_step(),
                                  125,
                                  138);
        break;
      case 6:
        check_solver_within_range((void)true,
                                  solver_control.last_step(),
                                  44,
                                  51);
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  deallog << "Eigenvalues:";
  for (unsigned int i = 0; i < v.size(); ++i)
    {
      deallog << ' ' << v[i];
      if (i != (v.size() - 1))
        deallog << ',';
    }
  deallog << std::endl << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(6);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    SolverControl control(500,
                          1e-12 /*1000*PETSC_MACHINE_EPSILON*/,
                          false,
                          false);


    const unsigned int size = 31;
    unsigned int       dim  = (size - 1);

    const unsigned n_eigenvalues = 4;

    deallog << "Size " << size << " Unknowns " << dim << std::endl << std::endl;

    // Make matrix
    FD1DLaplaceMatrix           testproblem(size);
    PETScWrappers::SparseMatrix A(dim, dim, 3);
    testproblem.three_point(A);
    A.compress(VectorOperation::insert);

    std::vector<PETScWrappers::MPI::Vector> u(
      n_eigenvalues, PETScWrappers::MPI::Vector(MPI_COMM_WORLD, dim, dim));
    std::vector<PetscScalar> v(n_eigenvalues);

    {
      SLEPcWrappers::SolverKrylovSchur solver(control);
      check_solve(solver, 1, control, A, u, v);
    }

    {
      SLEPcWrappers::SolverArnoldi solver(control);
      check_solve(solver, 2, control, A, u, v);
    }

    {
      SLEPcWrappers::SolverLanczos solver(control);
      check_solve(solver, 3, control, A, u, v);
    }

    {
      SLEPcWrappers::SolverGeneralizedDavidson solver(control);
      check_solve(solver, 4, control, A, u, v);
    }

    {
      SLEPcWrappers::SolverGeneralizedDavidson::AdditionalData data(true);
      SLEPcWrappers::SolverGeneralizedDavidson                 solver(control,
                                                      PETSC_COMM_SELF,
                                                      data);
      check_solve(solver, 5, control, A, u, v);
    }

    // set extra settings for JD; Otherwise, at least on OSX,
    // the number of eigensolver iterations is different between debug and
    // release modes!
    PETScWrappers::set_option_value("-st_ksp_type", "cg");
    PETScWrappers::set_option_value("-st_pc_type", "jacobi");
    PETScWrappers::set_option_value("-st_ksp_max_it", "10");
    {
      SLEPcWrappers::SolverJacobiDavidson solver(control);
      check_solve(solver, 6, control, A, u, v);
    }
  }
}
