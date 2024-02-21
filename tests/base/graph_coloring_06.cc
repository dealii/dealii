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



// Output the sizes of colors


#include <deal.II/base/graph_coloring.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <set>
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
  deallog << "dim=" << dim << std::endl;

  // Create the Triangulation and the DoFHandler
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // Create an adapted mesh
  for (unsigned int l = 0; l < 11 - 2 * dim; ++l)
    {
      typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active();
      for (; cell < dof_handler.end(); ++cell)
        if (cell->center().distance(Point<dim>()) < cell->diameter())
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }
  dof_handler.distribute_dofs(fe);

  deallog << "Total number of cells = " << triangulation.n_active_cells()
          << std::endl;

  // Create the coloring
  std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator>>
    coloring = GraphColoring::make_graph_coloring(
      dof_handler.begin_active(),
      dof_handler.end(),
      std::function<std::vector<types::global_dof_index>(
        const typename DoFHandler<dim>::active_cell_iterator &)>(
        &get_conflict_indices_cfem<dim>));

  for (unsigned int color = 0; color < coloring.size(); ++color)
    deallog << coloring[color].size() << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);
  deallog << std::fixed;

  check<1>();
  check<2>();
  check<3>();

  return 0;
}
