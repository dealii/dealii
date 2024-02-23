// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check that graph coloring works on adapted mesh with a uniform polynomial
// order.


#include <deal.II/base/graph_coloring.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>

#include "../tests.h"

template <int dim>
std::vector<types::global_dof_index>
get_conflict_indices_cfem(
  const typename DoFHandler<dim>::active_cell_iterator &it)
{
  std::vector<types::global_dof_index> local_dof_indices(
    it->get_fe().dofs_per_cell);
  it->get_dof_indices(local_dof_indices);

  return local_dof_indices;
}

template <int dim>
void
check()
{
  // Create the Triangulation and the DoFHandler
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // Create an adapted mesh
  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  for (; cell < dof_handler.end(); ++cell)
    {
      if ((cell->center()[0] == 0.625) && (cell->center()[1] == 0.625))
        cell->set_refine_flag();
    }
  triangulation.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);

  // Create the coloring
  std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator>>
    coloring(GraphColoring::make_graph_coloring(
      dof_handler.begin_active(),
      dof_handler.end(),
      std::function<std::vector<types::global_dof_index>(
        const typename DoFHandler<dim>::active_cell_iterator &)>(
        &get_conflict_indices_cfem<dim>)));

  // Output the coloring
  for (unsigned int color = 0; color < coloring.size(); ++color)
    {
      deallog << "Color: " << color << std::endl;
      for (unsigned int i = 0; i < coloring[color].size(); ++i)
        for (unsigned int j = 0; j < dim; ++j)
          deallog << coloring[color][i]->center()[j] << ' ';
      deallog << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);
  deallog << std::fixed;

  check<2>();

  return 0;
}
