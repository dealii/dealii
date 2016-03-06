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

// Tests for block_diagonal_operator

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
  FESystem<dim> fe(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1, FE_Q<dim>(3), 1);
  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  std::vector<types::global_dof_index> dpc (3);
  DoFTools::count_dofs_per_component (dof_handler, dpc);

  BlockDynamicSparsityPattern dsp(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      dsp.block(i, j).reinit (dpc[i], dpc[j]);
  dsp.collect_sizes ();

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  BlockSparseMatrix<double> a (sparsity_pattern);

  // Come up with a simple structure:

  for (unsigned int i = 0; i < a.block(0, 0).m(); ++i)
    a.block(0, 0).set(i, i, 10.);
  for (unsigned int i = 0; i < a.block(1, 1).m(); ++i)
    a.block(1, 1).set(i, i, 5.);
  for (unsigned int i = 0; i < a.block(2, 2).m(); ++i)
    a.block(2, 2).set(i, i, 3.);


  auto op_a = block_operator(a);

  auto op_b0 = linear_operator(a.block(0, 0));
  auto op_b1 = linear_operator(a.block(1, 1));
  auto op_b2 = linear_operator(a.block(2, 2));

  std::array<decltype(op_b0), 3> temp {{op_b0, op_b1, op_b2}};
  auto op_b = block_diagonal_operator<3, BlockVector<double>>(temp);


  {

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
    PRINTME("(A-B).vmult", x);

    x = 0.;
    op_x.vmult_add(x, u);
    PRINTME("(A-B).vmult_add", x);

    op_x.Tvmult(x, u);
    PRINTME("(A-B).Tvmult", x);

    x = 0.;
    op_x.Tvmult_add(x, u);
    PRINTME("(A-B).Tvmult_add", x);


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

  }

  // And finally the other block_diagonal_operator variant:

  std::array<decltype(op_b0), 5> temp2 {{op_b0, op_b0, op_b0, op_b0, op_b0}};
  auto op_c = block_diagonal_operator<5, BlockVector<double>>(temp2);

  auto op_d = block_diagonal_operator<5, BlockVector<double>>(op_b0);

  {
    BlockVector<double> u;
    op_c.reinit_domain_vector(u, false);
    for (unsigned int i = 0; i < u.size(); ++i)
      {
        u[i] = (double)(i+1);
      }
    PRINTME("u", u);

    BlockVector<double> x;
    op_c.reinit_range_vector(x, false);

    auto op_x = op_c - op_d;

    op_x.vmult(x, u);
    PRINTME("(C-D) vmult", x);

    x = 0.;
    op_x.vmult_add(x, u);
    PRINTME("(C-D) vmult_add", x);

    op_x.Tvmult(x, u);
    PRINTME("(C-D) Tvmult", x);

    x = 0.;
    op_x.Tvmult_add(x, u);
    PRINTME("(C-D) Tvmult_add", x);
  }

}
