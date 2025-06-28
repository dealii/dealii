// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the cell matrices generated in FETools and the local renumbering vector.

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <iostream>

#include "../tests.h"

void
logvec(const std::vector<types::global_dof_index> &v,
       const std::vector<types::global_dof_index> &w)
{
  deallog << '[';
  for (unsigned int i = 0; i < w.size(); ++i)
    deallog << ' ' << w[i];
  deallog << " ]" << std::endl << '[';
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << " ]" << std::endl;
}


template <int dim>
void
test_renumbering(const FiniteElement<dim> &fe)
{
  std::vector<types::global_dof_index> v(fe.dofs_per_cell);
  std::vector<types::global_dof_index> w(fe.n_blocks());

  deallog << fe.get_name() << std::endl;

  FETools::compute_block_renumbering(fe, v, w, false);
  logvec(v, w);
  FETools::compute_block_renumbering(fe, v, w, true);
  logvec(v, w);
}


template <int dim>
void
test_renumbering()
{
  deallog.push("Renumber");

  FE_Q<dim>             q1(1);
  FE_Q<dim>             q3(3);
  FE_DGQ<dim>           dg1(1);
  FE_RaviartThomas<dim> rt1(1);
  FE_Nedelec<dim>       n1(1);
  test_renumbering(q1);
  test_renumbering(q3);
  test_renumbering(dg1);
  test_renumbering(rt1);
  test_renumbering(n1);

  FESystem<dim> q1_3(q1, 3);
  test_renumbering(q1_3);
  FESystem<dim> q3_3(q3, 3);
  test_renumbering(q3_3);
  FESystem<dim> q3_2_q1_3(q3, 2, q1, 3);
  test_renumbering(q3_2_q1_3);
  test_renumbering(FESystem<dim>(rt1, 1, dg1, 1));
  test_renumbering(FESystem<dim>(rt1, 1, q1, 1));
  test_renumbering(FESystem<dim>(rt1, 2, q1, 2));
  test_renumbering(FESystem<dim>(n1, 1, q1, 1));

  deallog.pop();
}


int
main()
{
  initlog();
  //  test_renumbering<1>();
  test_renumbering<2>();
  test_renumbering<3>();
}
