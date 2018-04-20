// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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



// a un-hp-ified version of hp/step-2

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/dofs/dof_renumbering.h>



std::ofstream logfile("output");


void make_grid (Triangulation<2> &triangulation)
{
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
                              center, inner_radius, outer_radius,
                              10);
  //triangulation.reset_all_manifolds();
  GridTools::copy_boundary_to_manifold_id(triangulation);

  static const SphericalManifold<2> boundary_description(center);
  triangulation.set_manifold (0, boundary_description);

  for (unsigned int step=0; step<5; ++step)
    {
      Triangulation<2>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();

      for (; cell!=endc; ++cell)
        for (unsigned int vertex=0;
             vertex < GeometryInfo<2>::vertices_per_cell;
             ++vertex)
          {
            const double distance_from_center
              = center.distance (cell->vertex(vertex));

            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag ();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement ();
    }
}


void distribute_dofs (DoFHandler<2> &dof_handler)
{
  static const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs (finite_element);

  SparsityPattern sparsity_pattern (dof_handler.n_dofs(),
                                    dof_handler.n_dofs(),
                                    20);

  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();

  sparsity_pattern.print_gnuplot (deallog.get_file_stream());
}



void renumber_dofs (DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee (dof_handler);
  SparsityPattern sparsity_pattern (dof_handler.n_dofs(),
                                    dof_handler.n_dofs());

  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();

  sparsity_pattern.print_gnuplot (deallog.get_file_stream());
}





int main ()
{
  deallog << std::setprecision(2);

  deallog.attach(logfile);

  Triangulation<2> triangulation;
  make_grid (triangulation);

  DoFHandler<2> dof_handler (triangulation);

  distribute_dofs (dof_handler);
  renumber_dofs (dof_handler);
}
