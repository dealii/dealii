// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test TransferPrebuilt with system elements

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <algorithm>

#include "../tests.h"



template <typename Transfer>
void
make_matrix(const Transfer     &transfer,
            const unsigned int  high_level,
            FullMatrix<double> &matrix)
{
  Vector<double> src(matrix.n());
  Vector<double> dst(matrix.m());
  for (unsigned int i = 0; i < src.size(); ++i)
    {
      src    = 0;
      src(i) = 1;
      transfer.prolongate(high_level, dst, src);
      for (unsigned int j = 0; j < dst.size(); ++j)
        matrix(j, i) = dst(j);
    }
}



void
print_matrix(const FullMatrix<double> &m)
{
  for (unsigned int i = 0; i < m.m(); ++i)
    {
      for (unsigned int j = 0; j < m.n(); ++j)
        deallog << m(i, j) << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  DoFHandler<dim> mg_dof_handler(tr);
  mg_dof_handler.distribute_dofs(fe);
  mg_dof_handler.distribute_mg_dofs();

  DoFRenumbering::component_wise(mg_dof_handler);
  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    DoFRenumbering::component_wise(mg_dof_handler, level);

  MGTransferPrebuilt<Vector<double>> transfer;
  transfer.build(mg_dof_handler);

  FullMatrix<double> prolong_0_1(mg_dof_handler.n_dofs(1),
                                 mg_dof_handler.n_dofs(0));
  FullMatrix<double> prolong_1_2(mg_dof_handler.n_dofs(2),
                                 mg_dof_handler.n_dofs(1));

  deallog << "Level 0->1" << std::endl;
  make_matrix(transfer, 1, prolong_0_1);
  print_matrix(prolong_0_1);

  deallog << std::endl;

  deallog << "Level 1->2" << std::endl;
  make_matrix(transfer, 2, prolong_1_2);
  print_matrix(prolong_1_2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  // TODO: do in 1d
  check(FESystem<2>(FE_Q<2>(1), 2));
}
