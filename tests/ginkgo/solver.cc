// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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
#include <deal.II/lac/ginkgo_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <ginkgo/core/preconditioner/jacobi.hpp>

#include <iostream>
#include <typeinfo>

#include "../tests.h"

#include "../testmatrix.h"

int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);

  // std::shared_ptr<gko::Executor> exec = gko::CudaExecutor::create(0,
  //                                                               gko::OmpExecutor::create());
  // std::shared_ptr<gko::Executor> exec = gko::OmpExecutor::create();
  auto                           executor = "reference";
  std::shared_ptr<gko::Executor> exec     = gko::ReferenceExecutor::create();

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
    dealii::GinkgoWrappers::Csr<double> A(exec, dim, dim);
    testproblem.five_point(A);
    A.compress();

    GinkgoWrappers::Vector<double> f(exec, dim);
    f = 1.;
    GinkgoWrappers::Vector<double> u(exec, dim);
    u = 0.;

    GinkgoWrappers::SolverCG<>       cg_solver(exec, control);
    GinkgoWrappers::SolverBicgstab<> bicgstab_solver(exec, control);
    GinkgoWrappers::SolverCGS<>      cgs_solver(exec, control);
    GinkgoWrappers::SolverFCG<>      fcg_solver(exec, control);
    GinkgoWrappers::SolverGMRES<>    gmres_solver(exec, control);
    GinkgoWrappers::SolverIR<>       ir_solver_cg(exec, control);
    GinkgoWrappers::SolverCG<> cg_solver_with_jacobi_precond(exec, control);

    GinkgoWrappers::PreconditionBase<double, int> inner_cg(
      A,
      gko::solver::Cg<>::build()
        .with_criteria(gko::stop::Iteration::build().with_max_iters(45u).on(
                         exec),
                       gko::stop::ResidualNormReduction<>::build()
                         .with_reduction_factor(1e-5)
                         .on(exec))
        .on(exec));
    GinkgoWrappers::PreconditionJacobi<double, int> jacobi(A);

    check_solver_within_range(cg_solver.solve(A, u, f),
                              control.last_step(),
                              35,
                              39);
    u = 0.;
    check_solver_within_range(bicgstab_solver.solve(A, u, f),
                              control.last_step(),
                              26,
                              65);
    u = 0.;
    check_solver_within_range(cgs_solver.solve(A, u, f),
                              control.last_step(),
                              36,
                              79);
    u = 0.;
    check_solver_within_range(fcg_solver.solve(A, u, f),
                              control.last_step(),
                              33,
                              39);
    u = 0.;
    check_solver_within_range(gmres_solver.solve(A, u, f),
                              control.last_step(),
                              20,
                              49);
    u = 0.;
    check_solver_within_range(ir_solver_cg.solve(A, u, f, inner_cg),
                              control.last_step(),
                              0,
                              2);
    u = 0.;
    check_solver_within_range(
      cg_solver_with_jacobi_precond.solve(A, u, f, jacobi),
      control.last_step(),
      29,
      33);
  }
}
