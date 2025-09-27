// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// a un-hp-ified version of hp/step-2

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"



void
make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  // triangulation.reset_all_manifolds();
  GridTools::copy_boundary_to_manifold_id(triangulation);

  static const SphericalManifold<2> boundary_description(center);
  triangulation.set_manifold(0, boundary_description);

  for (unsigned int step = 0; step < 5; ++step)
    {
      Triangulation<2>::active_cell_iterator cell =
                                               triangulation.begin_active(),
                                             endc = triangulation.end();

      for (; cell != endc; ++cell)
        for (const unsigned int vertex : GeometryInfo<2>::vertex_indices())
          {
            const double distance_from_center =
              center.distance(cell->vertex(vertex));

            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement();
    }
}


void
distribute_dofs(DoFHandler<2> &dof_handler)
{
  static const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  SparsityPattern sparsity_pattern(dof_handler.n_dofs(),
                                   dof_handler.n_dofs(),
                                   20);

  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  sparsity_pattern.print_gnuplot(deallog.get_file_stream());
}



void
renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);
  SparsityPattern sparsity_pattern(dof_handler.n_dofs(), dof_handler.n_dofs());

  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  sparsity_pattern.print_gnuplot(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog << std::setprecision(2);


  Triangulation<2> triangulation;
  make_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
