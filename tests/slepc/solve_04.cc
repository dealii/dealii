// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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


// test the SLEPc solvers with HEP
// the unit tests is a mirror of
// SLEPc-3.6.1/src/eps/examples/tests/test4.c


#include "../tests.h"
#include "../lac/testmatrix.h"
#include "testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_compatibility.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <typeinfo>

template<typename SolverType, typename MatrixType, typename VectorType>
void
check_solve( SolverType &solver,
             const SolverControl &solver_control,
             const MatrixType &A,
             std::vector<VectorType> &u, std::vector<PetscScalar > &v)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  try
    {
      solver.set_problem_type(EPS_HEP);
      // reset vectors and set them as initial space
      // to avoid dependency on random numbers:
      for (unsigned int i = 0; i < u.size(); i++)
        for (unsigned int j=0; j<u[i].size(); ++j)
          u[i][j] = static_cast<double>(Testing::rand())/static_cast<double>(RAND_MAX);

      solver.set_initial_space(u);
      // solve
      solver.solve(A,v,u,v.size());
    }
  catch (dealii::SolverControl::NoConvergence &e)
    {
      deallog << "Exception: " << e.get_exc_name() << std::endl;
    }

  deallog << "Solver stopped after " << solver_control.last_step()
          << " iterations" << std::endl;

  deallog << "Eigenvalues:";
  for (unsigned int i = 0; i < v.size(); i++)
    {
      deallog << " "<<v[i];
      if ( i != (v.size()-1) )
        deallog <<",";
    }
  deallog << std::endl << std::endl;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(6);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  {
    SolverControl control(500,1e-12 /*1000*PETSC_MACHINE_EPSILON*/,false,false);


    const unsigned int size = 31;
    unsigned int dim = (size-1);

    const unsigned n_eigenvalues = 4;

    deallog << "Size " << size << " Unknowns " << dim << std::endl << std::endl;

    // Make matrix
    FD1DLaplaceMatrix testproblem(size);
    PETScWrappers::SparseMatrix  A(dim, dim, 3);
    testproblem.three_point(A);
    A.compress (VectorOperation::insert);

    std::vector<PETScWrappers::Vector>  u(n_eigenvalues,
                                          PETScWrappers::Vector(dim));
    std::vector<PetscScalar> v(n_eigenvalues);

    {
      SLEPcWrappers::SolverKrylovSchur solver(control);
      check_solve (solver, control, A,u,v);
    }

    {
      SLEPcWrappers::SolverArnoldi solver(control);
      check_solve (solver, control, A,u,v);
    }

    {
      SLEPcWrappers::SolverLanczos solver(control);
      check_solve (solver, control, A,u,v);
    }

    {
      SLEPcWrappers::SolverGeneralizedDavidson solver(control);
      check_solve (solver, control, A,u,v);
    }

    {
      SLEPcWrappers::SolverGeneralizedDavidson::AdditionalData data(true);
      SLEPcWrappers::SolverGeneralizedDavidson solver(control,PETSC_COMM_SELF,data);
      check_solve (solver, control, A,u,v);
    }

    // set extra settings for JD; Otherwise, at least on OSX,
    // the number of eigensolver iterations is different between debug and release modes!
    PETScWrappers::set_option_value("-st_ksp_type","cg");
    PETScWrappers::set_option_value("-st_pc_type", "jacobi");
    PETScWrappers::set_option_value("-st_ksp_max_it", "10");
    {
      SLEPcWrappers::SolverJacobiDavidson solver(control);
      check_solve (solver, control, A,u,v);
    }

  }

}

