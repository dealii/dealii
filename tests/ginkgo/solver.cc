// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include "../tests.h"

#include "../testmatrix.h"

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  {
    SolverControl control(200, 1e-3);

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
    auto                               executor = "reference";
    std::shared_ptr<gko::Executor>     exec = gko::ReferenceExecutor::create();
    std::shared_ptr<gko::LinOpFactory> jacobi =
      gko::preconditioner::Jacobi<>::build().on(exec);
    std::shared_ptr<gko::LinOpFactory> inner_cg =
      gko::solver::Cg<>::build()
        .with_criteria(gko::stop::Iteration::build().with_max_iters(45u).on(
                         exec),
                       gko::stop::ResidualNormReduction<>::build()
                         .with_reduction_factor(1e-5)
                         .on(exec))
        .on(exec);
    GinkgoWrappers::SolverCG<>       cg_solver(control, executor);
    GinkgoWrappers::SolverBicgstab<> bicgstab_solver(control, executor);
    GinkgoWrappers::SolverCGS<>      cgs_solver(control, executor);
    GinkgoWrappers::SolverFCG<>      fcg_solver(control, executor);
    GinkgoWrappers::SolverGMRES<>    gmres_solver(control, executor);
    GinkgoWrappers::SolverIR<>       ir_solver_cg(control, executor, inner_cg);
    GinkgoWrappers::SolverCG<>       cg_solver_with_jacobi_precond(control,
                                                             executor,
                                                             jacobi);

    check_solver_within_range(cg_solver.solve(A, u, f),
                              control.last_step(),
                              35,
                              39);
    u = 0.;
    check_solver_within_range(bicgstab_solver.solve(A, u, f),
                              control.last_step(),
                              53,
                              65);
    u = 0.;
    check_solver_within_range(cgs_solver.solve(A, u, f),
                              control.last_step(),
                              72,
                              79);
    u = 0.;
    check_solver_within_range(fcg_solver.solve(A, u, f),
                              control.last_step(),
                              33,
                              39);
    u = 0.;
    check_solver_within_range(gmres_solver.solve(A, u, f),
                              control.last_step(),
                              23,
                              29);
    u = 0.;
    check_solver_within_range(ir_solver_cg.solve(A, u, f),
                              control.last_step(),
                              0,
                              2);
    u = 0.;
    check_solver_within_range(cg_solver_with_jacobi_precond.solve(A, u, f),
                              control.last_step(),
                              29,
                              33);
  }
}
