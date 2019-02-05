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

// test the Ginkgo CG solver

#include <deal.II/lac/ginkgo_solver.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <iostream>
#include <typeinfo>

#include "../testmatrix.h"
#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  {
    SolverControl control(100, 1e-3);

    const unsigned int size = 32;
    unsigned int       dim  = (size - 1) * (size - 1);

    deallog << "Size " << size << " Unknowns " << dim << std::endl;

    // Make matrix
    FDMatrix        testproblem(size, size);
    SparsityPattern structure(dim, dim, 5);
    testproblem.five_point_structure(structure);
    structure.compress();
    SparseMatrix<double> A(structure);
    testproblem.five_point(A);

    Vector<double> f(dim);
    f = 1.;
    Vector<double> u(dim);
    u = 0.;

    // std::shared_ptr<gko::Executor> exec = gko::CudaExecutor::create(0,
    //                                                               gko::OmpExecutor::create());
    // std::shared_ptr<gko::Executor> exec = gko::OmpExecutor::create();
    auto executor = "reference";
    std::shared_ptr<gko::Executor> exec = gko::ReferenceExecutor::create();
    std::shared_ptr<gko::LinOpFactory> jacobi = gko::preconditioner::Jacobi<>::build().on(exec);
    GinkgoWrappers::SolverCG<>     solver(control, executor);
    GinkgoWrappers::SolverCG<>     solver_with_jacobi_precond(control, executor, jacobi);

    check_solver_within_range(solver.solve(A, u, f),
                              control.last_step(),
                              35,
                              39);
    u = 0.;
    check_solver_within_range(solver_with_jacobi_precond.solve(A, u, f),
                              control.last_step(),
                              29,
                              33);
  }
}
