// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests for the LinearOperator template with
//   dealii::Vector<double>
//   dealii::SparseMatrix<double>
//   dealii::FullMatrix<double>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  static const int dim = 2;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  MappingQ<dim>   mapping_q1(1);
  FE_Q<dim>       q1(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(q1);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  dsp.compress();
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  SparseMatrix<double> a(sparsity_pattern);
  SparseMatrix<double> b(sparsity_pattern);

  QGauss<dim> quadrature(4);
  MatrixCreator::create_laplace_matrix(mapping_q1, dof_handler, quadrature, a);
  MatrixCreator::create_mass_matrix(mapping_q1, dof_handler, quadrature, b);


  // Constructors and assignment:

  auto op_a = linear_operator(a);
  auto op_b = linear_operator(b);

  {
    LinearOperator<dealii::Vector<double>, dealii::Vector<double>> op_x(a);
    op_a = a;
    op_b = b;
  }

  // vmult:

  Vector<double> u;
  op_a.reinit_domain_vector(u, true);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i + 1);
    }

  deallog << "u: " << u << std::endl;

  Vector<double> v;
  op_a.reinit_domain_vector(v, false);
  Vector<double> w;
  op_a.reinit_domain_vector(w, false);
  Vector<double> x;
  op_a.reinit_domain_vector(x, false);

  op_a.vmult(v, u);
  deallog << "Au: " << v << std::endl;

  op_b.vmult(w, u);
  deallog << "Bu: " << w << std::endl;

  // operator+, operator-, operator+=, operator-=:

  x = v;
  x += w;
  deallog << "Au+Bu: " << x << std::endl;

  (op_a + op_b).vmult(x, u);
  deallog << "(A+B)u: " << x << std::endl;

  auto op_x = op_a;
  op_x += op_b;
  op_x.vmult(x, u);
  deallog << "(A+=B)u: " << x << std::endl;

  x = v;
  x -= w;
  deallog << "Au-Bu: " << x << std::endl;

  (op_a - op_b).vmult(x, u);
  deallog << "(A-B)u: " << x << std::endl;

  op_x = op_a;
  op_x -= op_b;
  op_x.vmult(x, u);
  deallog << "(A-=B)u: " << x << std::endl;

  // operator*, operator*=

  op_b.vmult(v, u);
  op_a.vmult(w, v);
  deallog << "(A(Bu)): " << w << std::endl;

  (op_a * op_b).vmult(x, u);
  deallog << "(A*B)u: " << x << std::endl;

  op_x = op_a;
  op_x *= op_b;
  op_x.vmult(x, u);
  deallog << "(A*=B)u: " << x << std::endl;

  op_x *= 4.;
  op_x.vmult(x, u);
  deallog << "(A*=B*=4.)u: " << x << std::endl;

  // solver interface:

  SolverControl solver_control(1000, 1e-12);
  SolverCG<>    solver(solver_control);

  deallog.depth_file(0);
  solver.solve(op_b, v, u, PreconditionIdentity());
  deallog.depth_file(3);
  deallog << "solve(B, v, u): " << v << std::endl;

  deallog.depth_file(0);
  inverse_operator(op_b, solver, PreconditionIdentity()).vmult(v, u);
  deallog.depth_file(3);
  deallog << "inverse_operator(B)u: " << v << std::endl;

  deallog.depth_file(0);
  op_b.vmult(w, v);
  deallog.depth_file(3);
  deallog << "B(inverse_operator(B)u): " << w << std::endl;

  deallog.depth_file(0);
  (op_b * inverse_operator(op_b, solver, PreconditionIdentity())).vmult(w, u);
  deallog.depth_file(3);
  deallog << "(B*inverse_operator(B))u: " << w << std::endl;

  deallog.depth_file(0);
  (inverse_operator(op_b, solver, PreconditionIdentity()) * op_b).vmult(w, u);
  deallog.depth_file(3);
  deallog << "(inverse_operator(B)*B)u: " << w << std::endl;

  SolverControl inner_solver_control(1000, 1e-12);
  SolverCG<>    inner_solver(solver_control);

  deallog.depth_file(0);
  solver.solve(inverse_operator(op_b, inner_solver, PreconditionIdentity()),
               v,
               u,
               PreconditionIdentity());
  deallog.depth_file(3);
  deallog << "solve(inverse_operator(B), v, u) == Bu: " << v << std::endl;

  deallog.depth_file(0);
  solver.solve(op_b + 0.005 * op_a, v, u, PreconditionIdentity());
  deallog.depth_file(3);
  deallog << "solve(B+0.5*0.01*A, v, u): " << v << std::endl;
}
