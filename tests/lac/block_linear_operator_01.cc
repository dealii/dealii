// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

// Tests for block_operator

#include "../tests.h"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#define PRINTME(name, var) \
  deallog << "Block vector: " name << ":" << std::endl; \
  for (unsigned int i = 0; i < var.n_blocks(); ++i) \
    deallog << "[block " << i << " ]  " << var.block(i);


using namespace dealii;

int main()
{
  initlog();
  deallog << std::setprecision(10);

  static const int dim = 2;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global(2);

  MappingQGeneric<dim> mapping_q1(1);
  FESystem<dim> fe(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  std::vector<types::global_dof_index> dofs_per_component (2);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
  const unsigned int n_u = dofs_per_component[0],
                     n_p = dofs_per_component[1];

  BlockDynamicSparsityPattern dsp(2, 2);
  dsp.block(0, 0).reinit (n_u, n_u);
  dsp.block(1, 0).reinit (n_p, n_u);
  dsp.block(0, 1).reinit (n_u, n_p);
  dsp.block(1, 1).reinit (n_p, n_p);

  for (unsigned int i = 0; i < n_u; ++i)
    {
      for (unsigned int j = 0; j < n_u; ++j)
        dsp.block(0, 0).add(i, j);
      for (unsigned int j = 0; j < n_p; ++j)
        {
          dsp.block(0, 1).add(i, j);
          dsp.block(1, 0).add(j, i);
        }
    }
  for (unsigned int i = 0; i < n_p; ++i)
    for (unsigned int j = 0; j < n_p; ++j)
      dsp.block(1, 1).add(i, j);
  dsp.collect_sizes ();

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  BlockSparseMatrix<double> a (sparsity_pattern);

  // Come up with some simple, but non-trivial structure:

  for (unsigned int i = 0; i < a.block(0, 0).n(); ++i)
    a.block(0, 0).set(i, i, 10.);

  for (unsigned int i = 0; i < a.block(1, 1).n(); ++i)
    a.block(1, 1).set(i, i, 5.);

  for (unsigned int i = 0; i < a.block(1, 0).m(); ++i)
    a.block(0, 1).set(i, i, 3.);

  for (unsigned int i = 0; i < a.block(0, 1).n(); ++i)
    a.block(1, 0).set(i, i, 1.);


  auto op_a = linear_operator<BlockVector<double>>(a);

  auto op_b00 = linear_operator(a.block(0, 0));
  auto op_b01 = linear_operator(a.block(0, 1));
  auto op_b10 = linear_operator(a.block(1, 0));
  auto op_b11 = linear_operator(a.block(1, 1));

  std::array<std::array<decltype(op_b00), 2>, 2> temp
  {
    {op_b00, op_b01, op_b10, op_b11}
  };
  auto op_b = block_operator<2, 2, BlockVector<double>>(temp);

  {
    // also test copy and reference assignment to a LinearOperator
    LinearOperator<BlockVector<double>> &op_x1 = op_b;
    LinearOperator<BlockVector<double>>  op_x2 = op_b;
    LinearOperator<BlockVector<double>> op_x3(op_b);
  }

  // vmult:

  BlockVector<double> u;
  op_a.reinit_domain_vector(u, false);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i+1);
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

  x = v;
  x -= w;
  PRINTME("Au-Bu", x);

  // Test that both objects give the same results:

  auto op_x = op_a - op_b;

  op_x.vmult(x, u);
  PRINTME("vmult", x);

  x = 0.;
  op_x.vmult_add(x, u);
  PRINTME("vmult_add", x);

  op_x.Tvmult(x, u);
  PRINTME("Tvmult", x);

  x = 0.;
  op_x.Tvmult_add(x, u);
  PRINTME("Tvmult_add", x);


  // Test vector reinitalization:

  op_x = op_b * op_b * op_b;
  op_x.vmult(x, u);
  PRINTME("(B*B*B) vmult", x);

  x = 0.;
  op_x.vmult_add(x, u);
  PRINTME("(B*B*B) vmult_add", x);

  op_x.Tvmult(x, u);
  PRINTME("(B*B*B) Tvmult", x);

  x = 0.;
  op_x.Tvmult_add(x, u);
  PRINTME("(B*B*B) Tvmult_add", x);

  // And finally complicated block structures:

  std::array<std::array<decltype(op_b00), 3>, 3> temp2
  {
    {op_b00, op_b01, op_b00, op_b10, op_b11, op_b10, op_b10, op_b11, op_b10}
  };
  auto op_upp_x_upu = block_operator<3, 3, BlockVector<double>>(temp2);

  op_upp_x_upu.reinit_domain_vector(u, false);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i+1);
    }
  PRINTME("u", u);

  op_upp_x_upu.reinit_range_vector(v, false);
  op_upp_x_upu.vmult(v, u);
  PRINTME("v", v);

  v = 0.;
  op_upp_x_upu.vmult_add(v, u);
  PRINTME("v", v);

  std::array<std::array<decltype(op_b01), 1>, 3> temp3
  {
    {op_b01, op_b11, op_b11}
  };
  auto op_upp_x_p = block_operator<3, 1, BlockVector<double>>(temp3);

  std::array<std::array<decltype(op_b01), 3>, 1> temp4
  {
    {op_b00, op_b01, op_b00}
  };
  auto op_u_x_upu = block_operator<1, 3, BlockVector<double>>(temp4);

  auto op_long = op_u_x_upu * transpose_operator(op_upp_x_upu) * op_upp_x_p;

  op_long.reinit_domain_vector(u, false);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i+1);
    }
  PRINTME("u", u);

  op_long.reinit_range_vector(v, false);
  op_long.vmult(v, u);
  PRINTME("v", v);

  v = 0.;
  op_long.vmult_add(v, u);
  PRINTME("v", v);

}



