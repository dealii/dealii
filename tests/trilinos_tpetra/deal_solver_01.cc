// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the CG solver using the Trilinos matrix and vector classes


#include <deal.II/base/utilities.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_qmrs.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <typeinfo>

#include "../tests.h"

#include "../testmatrix.h"


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  {
    SolverControl control(200, 1.e-3);

    const unsigned int size = 32;
    unsigned int       dim  = (size - 1) * (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix               testproblem(size, size);
    DynamicSparsityPattern csp(dim, dim);
    testproblem.five_point_structure(csp);
    LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> A;
    A.reinit(csp);
    testproblem.five_point(A);

    LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> f;
    f.reinit(complete_index_set(dim), MPI_COMM_WORLD);
    LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> u;
    u.reinit(complete_index_set(dim), MPI_COMM_WORLD);
    f.trilinos_rcp()->putScalar(1.0);
    A.compress(VectorOperation::insert);
    f.compress(VectorOperation::insert);
    u.compress(VectorOperation::insert);

    GrowingVectorMemory<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
      mem;
    SolverCG<
      LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>>
                         solver(control, mem);
    PreconditionIdentity preconditioner;

    deallog << "Solver type: " << typeid(solver).name() << std::endl;

    check_solver_within_range(solver.solve(A, u, f, preconditioner),
                              control.last_step(),
                              42,
                              44);
  }
}
