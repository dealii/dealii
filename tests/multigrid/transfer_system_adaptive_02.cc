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


// test like _01 but with boundary conditions
#include <deal.II/base/function.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

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
refine_mesh(Triangulation<dim> &triangulation)
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    {
      const Point<dim> p        = cell->center();
      bool             positive = p[0] > 0;
      if (positive)
        {
          cell->set_refine_flag();
          cell_refined = true;
        }
    }
  if (!cell_refined) // if no cell was selected for refinement, refine global
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);

  std::vector<unsigned int> subdivisions(dim, 1);
  subdivisions[0] = 2;

  const Point<dim> bottom_left =
    (dim == 2 ? Point<dim>(-1, -1) : Point<dim>(-1, -1, -1));
  const Point<dim> top_right =
    (dim == 2 ? Point<dim>(1, 1) : Point<dim>(1, 1, 1));
  GridGenerator::subdivided_hyper_rectangle(
    tr, subdivisions, bottom_left, top_right, true);
  refine_mesh(tr);

  DoFHandler<dim> mg_dof_handler(tr);
  mg_dof_handler.distribute_dofs(fe);
  mg_dof_handler.distribute_mg_dofs();

  deallog << "Global  dofs: " << mg_dof_handler.n_dofs() << std::endl;
  for (unsigned int l = 0; l < tr.n_levels(); ++l)
    {
      deallog << "Level " << l << " dofs:";
      deallog << ' ' << mg_dof_handler.n_dofs(l);
      deallog << std::endl;
    }

  DoFRenumbering::component_wise(mg_dof_handler);
  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    DoFRenumbering::component_wise(mg_dof_handler, level);

  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mg_dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(mg_dof_handler, {3});

  MGTransferPrebuilt<Vector<double>> transfer(mg_constrained_dofs);
  transfer.build(mg_dof_handler);

  FullMatrix<double> prolong_0_1(mg_dof_handler.n_dofs(1),
                                 mg_dof_handler.n_dofs(0));

  deallog << "Level 0->1" << std::endl;
  make_matrix(transfer, 1, prolong_0_1);
  print_matrix(prolong_0_1);
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  // TODO: do in 1d
  check(FESystem<2>(FE_Q<2>(1), 2));
}
