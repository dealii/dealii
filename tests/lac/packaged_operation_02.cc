// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test PackagedOperation for
//   dealii::SparseMatrix<double> with LinearOperator

#include <deal.II/lac/packaged_operation.h>

#include "../tests.h"

// and a _lot_ of stuff to create a linera oprator
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



void
test_applies(std::string                              description,
             const PackagedOperation<Vector<double>> &expr)
{
  // test apply
  Vector<double> tmp = expr;
  deallog << description << ": " << tmp << std::endl;

  // test apply_add
  for (auto &i : tmp)
    i = 100.;
  expr.apply_add(tmp);
  deallog << "100. * 1_n + " << description << ": " << tmp << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  static const int dim = 2;

  // Create mass marix M, and an iterative inverse MInv:

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  MappingQGeneric<dim> mapping_q1(1);
  FE_Q<dim>            q1(1);
  DoFHandler<dim>      dof_handler(triangulation);
  dof_handler.distribute_dofs(q1);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  dsp.compress();
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  SparseMatrix<double> m(sparsity_pattern);

  QGauss<dim> quadrature(4);
  MatrixCreator::create_mass_matrix(mapping_q1, dof_handler, quadrature, m);

  auto M = linear_operator(m);

  SolverControl solver_control(1000, 1e-12);
  SolverCG<>    solver(solver_control);

  auto MInv = inverse_operator(M, solver, PreconditionIdentity());

  // Tests:

  Vector<double> u;
  M.reinit_domain_vector(u, true);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i + 1);
    }

  deallog << "u: " << u << std::endl;

  deallog.depth_file(0);
  Vector<double> b = MInv * u;
  deallog.depth_file(3);
  deallog << "b: " << b << std::endl;

  test_applies("M * b", M * b);
  test_applies("b * M", b * M);

  auto expr = b;
  test_applies("M * b", M * expr);
  test_applies("b * M", expr * M);
}
