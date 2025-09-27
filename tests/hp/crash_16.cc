// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// an extract of hp_constraints_q_system_x_01. something still goes wrong,
// akin to what happens in crash_15. (it turned out to be bogus index
// computations.)

char logname[] = "output";


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  hp::FECollection<dim> fe;
  for (unsigned int i = 1; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      fe.push_back(FESystem<dim>(FE_Q<dim>(i), 1, FE_DGQ<dim>(j), 1));

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  DoFHandler<dim> dof_handler(triangulation);

  // distribute fe_indices randomly
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    cell->set_active_fe_index(Testing::rand() % fe.size());
  dof_handler.distribute_dofs(fe);

  // loop over all lines and make sure that
  // all the DoF indices on these lines are
  // identical
  std::vector<types::global_dof_index> indices_1;
  std::vector<types::global_dof_index> indices_2;

  std::set<unsigned int> line_already_treated;

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
      if (line_already_treated.find(cell->line(l)->index()) ==
          line_already_treated.end())
        // line not yet treated
        {
          const typename DoFHandler<dim>::active_line_iterator line =
            cell->line(l);
          deallog << "line=" << line << std::endl;
          line_already_treated.insert(line->index());

          for (unsigned int f = 0; f < line->n_active_fe_indices(); ++f)
            {
              indices_1.resize(
                fe[line->nth_active_fe_index(f)].dofs_per_line +
                2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
              line->get_dof_indices(indices_1, line->nth_active_fe_index(f));

              deallog << "  fe index=" << line->nth_active_fe_index(f)
                      << ", indices=";
              for (unsigned int i = 0; i < indices_1.size(); ++i)
                deallog << indices_1[i] << ' ';

              deallog << std::endl;
            }

          for (unsigned int f = 0; f < line->n_active_fe_indices(); ++f)
            {
              indices_1.resize(
                fe[line->nth_active_fe_index(f)].dofs_per_line +
                2 * fe[line->nth_active_fe_index(f)].dofs_per_vertex);
              line->get_dof_indices(indices_1, line->nth_active_fe_index(f));
              for (unsigned int g = f + 1; g < line->n_active_fe_indices(); ++g)
                if (fe[line->nth_active_fe_index(f)].dofs_per_line ==
                    fe[line->nth_active_fe_index(g)].dofs_per_line)
                  {
                    indices_2.resize(
                      fe[line->nth_active_fe_index(g)].dofs_per_line +
                      2 * fe[line->nth_active_fe_index(g)].dofs_per_vertex);
                    line->get_dof_indices(indices_2,
                                          line->nth_active_fe_index(g));
                    Assert(indices_1 == indices_2, ExcInternalError());
                  }
            }
        }
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  test<3>();
}
