// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

// Tests for the LinearOperator template with
//   dealii::BlockVector<double>
//   dealii::BlockSparseMatrix<double>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"

#define PRINTME(name, var)                                            \
  deallog << name << ": [block 0] " << var.block(0) << "  [block 1] " \
          << var.block(1) << std::endl;



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  static const int dim = 2;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  MappingQGeneric<dim> mapping_q1(1);
  FESystem<dim>        fe(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1);
  DoFHandler<dim>      dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  const std::vector<types::global_dof_index> dofs_per_component =
    DoFTools::count_dofs_per_fe_component(dof_handler);
  const unsigned int n_u = dofs_per_component[0], n_p = dofs_per_component[1];

  BlockDynamicSparsityPattern dsp(2, 2);
  dsp.block(0, 0).reinit(n_u, n_u);
  dsp.block(1, 0).reinit(n_p, n_u);
  dsp.block(0, 1).reinit(n_u, n_p);
  dsp.block(1, 1).reinit(n_p, n_p);
  dsp.collect_sizes();
  DoFTools::make_sparsity_pattern(dof_handler, dsp);

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  BlockSparseMatrix<double> a(sparsity_pattern);
  BlockSparseMatrix<double> b(sparsity_pattern);

  for (unsigned int i = 0; i < a.n(); ++i)
    {
      a.set(i, i, 1.);
      b.set(i, i, 5.);
    }

  // Constructors and assignment:

  auto op_a = linear_operator<BlockVector<double>>(a);
  auto op_b = linear_operator<BlockVector<double>>(b);

  {
    decltype(op_a) op_x(a);
    op_a = a;
    op_b = b;
  }

  // vmult:

  BlockVector<double> u;
  op_a.reinit_domain_vector(u, true);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i + 1);
    }

  PRINTME("u", u);

  BlockVector<double> v;
  op_a.reinit_domain_vector(v, false);
  BlockVector<double> w;
  op_a.reinit_domain_vector(w, false);
  BlockVector<double> x;
  op_a.reinit_domain_vector(x, false);

  op_a.vmult(v, u);
  PRINTME("Au", v);

  op_b.vmult(w, u);
  PRINTME("Bu", w);

  // operator+, operator-, operator+=, operator-=:

  x = v;
  x += w;
  PRINTME("Au+Bu", x);

  (op_a + op_b).vmult(x, u);
  PRINTME("(A+B)u", x);

  auto op_x = op_a;
  op_x += op_b;
  op_x.vmult(x, u);
  PRINTME("(A+=B)u", x);

  x = v;
  x -= w;
  PRINTME("Au-Bu", x);

  (op_a - op_b).vmult(x, u);
  PRINTME("(A-B)u", x);

  op_x = op_a;
  op_x -= op_b;
  op_x.vmult(x, u);
  PRINTME("(A-=B)u", x);

  // operator*, operator*=

  op_b.vmult(v, u);
  op_a.vmult(w, v);
  PRINTME("(A(Bu))", x);

  (op_a * op_b).vmult(x, u);
  PRINTME("(A*B)u", x);

  op_x = op_a;
  op_x *= op_b;
  op_x.vmult(x, u);
  PRINTME("(A*=B)u", x);

  // solver interface:

  SolverControl                    solver_control(1000, 1e-10);
  SolverGMRES<BlockVector<double>> solver(solver_control);

  deallog.depth_file(0);
  solver.solve(op_b, v, u, PreconditionIdentity());
  deallog.depth_file(3);
  PRINTME("solve(B, v, u)", v);

  deallog.depth_file(0);
  inverse_operator(op_b, solver, PreconditionIdentity()).vmult(v, u);
  deallog.depth_file(3);
  PRINTME("inverse_operator(B)u", v);

  deallog.depth_file(0);
  op_b.vmult(w, v);
  deallog.depth_file(3);
  PRINTME("B(inverse_operator(B)u)", w);

  deallog.depth_file(0);
  (op_b * inverse_operator(op_b, solver, PreconditionIdentity())).vmult(w, u);
  deallog.depth_file(3);
  PRINTME("(B*inverse_operator(B))u", w);

  deallog.depth_file(0);
  (inverse_operator(op_b, solver, PreconditionIdentity()) * op_b).vmult(w, u);
  deallog.depth_file(3);
  PRINTME("(inverse_operator(B)*B)u", w);

  SolverControl                 inner_solver_control(1000, 1e-12);
  SolverCG<BlockVector<double>> inner_solver(solver_control);

  deallog.depth_file(0);
  solver.solve(inverse_operator(op_b, inner_solver, PreconditionIdentity()),
               v,
               u,
               PreconditionIdentity());
  deallog.depth_file(3);
  PRINTME("solve(inverse_operator(B), v, u) == Bu", v);

  deallog.depth_file(0);
  solver.solve(op_b + 0.005 * op_a, v, u, PreconditionIdentity());
  deallog.depth_file(3);
  PRINTME("solve(B+0.5*0.01*A, v, u)", v);
}
