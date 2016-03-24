// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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


// test the BiCGStab solver using the PETSc matrix and vector classes



#include "../tests.h"
#include "../lac/testmatrix.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <typeinfo>

template<typename SolverType, typename MatrixType, typename VectorType, class PRECONDITION>
void
check_solve(SolverType          &solver,
            const SolverControl &solver_control,
            const MatrixType    &A,
            VectorType          &u,
            VectorType          &f,
            const PRECONDITION  &P)
{
  deallog << "Solver type: " << typeid(solver).name() << std::endl;

  u = 0.;
  f = 1.;
  try
    {
      check_solver_within_range(
        solver.solve(A,u,f,P),
        solver_control.last_step(), 48, 51);
    }
  catch (std::exception &e)
    {
      deallog << e.what() << std::endl;
      abort ();
    }
}


int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  {
    SolverControl control(100, 1.e-10);

    const unsigned int size = 32;
    unsigned int dim = (size-1)*(size-1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix testproblem(size, size);
    PETScWrappers::SparseMatrix  A(dim, dim, 5);
    testproblem.five_point(A);
    A.compress (VectorOperation::insert);

    PETScWrappers::Vector  f(dim);
    PETScWrappers::Vector  u(dim);
    f = 1.;

    GrowingVectorMemory<PETScWrappers::Vector> mem;
    SolverBicgstab<PETScWrappers::Vector> solver(control,mem);
    PreconditionIdentity preconditioner;
    check_solve (solver, control, A,u,f, preconditioner);
  }
  GrowingVectorMemory<PETScWrappers::Vector>::release_unused_memory ();

}
