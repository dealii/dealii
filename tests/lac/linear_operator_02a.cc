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

// Variant of linear_operator_02 that uses std::complex<double> as
// value_type: Tests for the LinearOperator template with
//   dealii::Vector<complex>
//   dealii::SparseMatrix<complex>
//   dealii::FullMatrix<complex>

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
#include <deal.II/lac/sparse_matrix.templates.h>

#include <deal.II/numerics/matrix_tools.h>

#include <complex>

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

  // use a real valued matrix and rely on the template overloads of vmult,
  // etc. to also work with std::complex<double>
  SparseMatrix<double> a(sparsity_pattern);
  SparseMatrix<double> b(sparsity_pattern);

  QGauss<dim> quadrature(4);
  MatrixCreator::create_laplace_matrix(mapping_q1, dof_handler, quadrature, a);
  MatrixCreator::create_mass_matrix(mapping_q1, dof_handler, quadrature, b);


  // Constructors and assignment:

  auto op_a = linear_operator<dealii::Vector<std::complex<double>>>(a);
  auto op_b = linear_operator<dealii::Vector<std::complex<double>>>(b);

  {
    LinearOperator<dealii::Vector<std::complex<double>>> op_x(a);
    op_a = a;
    op_b = b;
  }

  // vmult:

  Vector<std::complex<double>> u;
  op_a.reinit_domain_vector(u, true);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (std::complex<double>)(i + 1);
    }

  deallog << "u: " << u << std::endl;

  Vector<std::complex<double>> v;
  op_a.reinit_domain_vector(v, false);
  Vector<std::complex<double>> w;
  op_a.reinit_domain_vector(w, false);
  Vector<std::complex<double>> x;
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

  op_x *= std::complex<double>(4.);
  op_x.vmult(x, u);
  deallog << "(A*=B*=4.)u: " << x << std::endl;

  // solver interface:

  // TODO: Well, update as soon as SolverCG can work with complex data
  // types.
}
