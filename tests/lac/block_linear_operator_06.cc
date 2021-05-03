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

// Tests for block_operator for operations with identical source destinations

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

#include "../tests.h"

#define PRINTME(name, var)                                            \
  deallog << "Block vector: " name << ":" << std::endl;               \
  for (unsigned int i = 0; i < var.n_blocks(); ++i)                   \
    deallog << "[block " << i << " ]  " << var.block(i) << std::endl; \
  deallog << std::endl;



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
  FESystem<dim>        fe(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1);
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
  dsp.collect_sizes();

  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dsp);
  sparsity_pattern.compress();

  BlockSparseMatrix<double> a(sparsity_pattern);

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

  std::array<std::array<decltype(op_b00), 2>, 2> temp{
    {{{op_b00, op_b01}}, {{op_b10, op_b11}}}};
  auto op_b = block_operator<2, 2, BlockVector<double>>(temp);

  // vmult:

  BlockVector<double> u;
  op_a.reinit_domain_vector(u, false);
  for (unsigned int i = 0; i < u.size(); ++i)
    {
      u[i] = (double)(i + 1);
    }
  BlockVector<double> u_copy;
  op_a.reinit_domain_vector(u_copy, false);
  u_copy = u;

  BlockVector<double> result;
  op_a.reinit_domain_vector(result, false);

  PRINTME("u", u);
  PRINTME("u_copy", u);
  {
    op_a.vmult(u, u);
    PRINTME("A.vmult", u);
    result = u;
    u      = u_copy;

    op_b.vmult(u, u);
    PRINTME("B.vmult", u);
    result -= u;
    u = u_copy;

    AssertThrow(result.linfty_norm() < 1.e-10, ExcInternalError());
  }

  {
    op_a.Tvmult(u, u);
    PRINTME("A.Tvmult", u);
    result = u;
    u      = u_copy;

    op_b.Tvmult(u, u);
    PRINTME("B.Tvmult", u);
    result -= u;
    u = u_copy;

    AssertThrow(result.linfty_norm() < 1.e-10, ExcInternalError());
  }

  {
    op_a.vmult_add(u, u);
    u -= u_copy;
    PRINTME("A.vmult_add", u);
    result = u;
    u      = u_copy;

    op_b.vmult_add(u, u);
    u -= u_copy;
    PRINTME("B.vmult_add", u);
    result -= u;
    u = u_copy;

    AssertThrow(result.linfty_norm() < 1.e-10, ExcInternalError());
  }

  {
    op_a.Tvmult_add(u, u);
    PRINTME("A.Tvmult_add", u);
    result = u;
    u      = u_copy;

    op_b.Tvmult_add(u, u);
    PRINTME("B.Tvmult_add", u);
    result -= u;

    AssertThrow(result.linfty_norm() < 1.e-10, ExcInternalError());
  }
}
