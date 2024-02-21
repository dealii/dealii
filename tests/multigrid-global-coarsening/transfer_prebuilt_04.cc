// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Deadlock reported by Kronbichler (github
// https://github.com/dealii/dealii/issues/2051) with 3 processes
// in MGTransferPrebuilt


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"


template <int dim>
void
check()
{
  FE_Q<dim> fe(1);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tr, 3);
  for (unsigned int cycle = 0; cycle < (dim == 2 ? 10 : 7); ++cycle)
    {
      // adaptive refinement into a circle
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tr.begin_active();
           cell != tr.end();
           ++cell)
        if (cell->is_locally_owned() && cell->vertex(0).norm() < 1e-10)
          cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      deallog << "no. cells: " << tr.n_global_active_cells() << " on "
              << tr.n_global_levels() << " levels" << std::endl;

      DoFHandler<dim> mgdof(tr);
      mgdof.distribute_dofs(fe);
      mgdof.distribute_mg_dofs();

      MGConstrainedDoFs mg_constrained_dofs;
      mg_constrained_dofs.initialize(mgdof);
      mg_constrained_dofs.make_zero_boundary_constraints(mgdof, {0});

      unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      if (0)
        {
          std::ofstream grid_output(
            ("out" + Utilities::to_string(myid) + ".svg").c_str());
          GridOut           grid_out;
          GridOutFlags::Svg flags;
          flags.label_level_subdomain_id = true;
          flags.coloring = GridOutFlags::Svg::level_subdomain_id;
          flags.convert_level_number_to_height = true;
          grid_out.set_flags(flags);

          grid_out.write_svg(tr, grid_output);
        }
#ifdef DEAL_II_WITH_TRILINOS
      {
        // MGTransferPrebuilt internally uses Trilinos matrices, so only
        // create this if we have Trilinos
        MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>>
          transfer_ref(mg_constrained_dofs);
        transfer_ref.build(mgdof);
      }
#endif
      {
        // but the matrix free transfer will work without Trilinos
        MGTransferMF<dim, double> transfer_ref(mg_constrained_dofs);
        transfer_ref.build(mgdof);
      }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  mpi_initlog();

  check<2>();
}
