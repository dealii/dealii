// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that DofRenumbering::Cuthill_McKee works in parallel (by applying
// Cuthill-McKee individually on each processor's subdomain).


#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int nprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tr, -1.0, 1.0);
  tr.refine_global(8 - 2 * dim);

  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      if (cell->center().norm() < 0.3)
        {
          cell->set_refine_flag();
        }

  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(1);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();

  for (unsigned int level = 0; level <= tr.n_levels(); ++level)
    {
      types::global_dof_index              n_dofs{0};
      std::vector<types::global_dof_index> renumbering;
      if (level < tr.n_levels())
        {
          n_dofs      = dofh.n_dofs(level);
          renumbering = std::vector<types::global_dof_index>(
            dofh.locally_owned_mg_dofs(level).n_elements());
          DoFRenumbering::compute_Cuthill_McKee(
            renumbering,
            dofh,
            false,
            false,
            std::vector<types::global_dof_index>(),
            level);
        }
      else
        {
          n_dofs      = dofh.n_dofs();
          renumbering = std::vector<types::global_dof_index>(
            dofh.locally_owned_dofs().n_elements());
          DoFRenumbering::compute_Cuthill_McKee(renumbering, dofh);
        }

      // send everything to processor 0 for output
      std::vector<types::global_dof_index> complete_renumbering(n_dofs);
      std::copy(renumbering.begin(),
                renumbering.end(),
                complete_renumbering.begin());
      unsigned int          offset = renumbering.size();
      std::vector<IndexSet> dofs_per_proc;
      if (level < tr.n_levels())
        {
          dofs_per_proc =
            Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                       dofh.locally_owned_mg_dofs(level));
        }
      else
        {
          dofs_per_proc = Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                                     dofh.locally_owned_dofs());
        }

      for (unsigned int i = 1; i < nprocs; ++i)
        {
          if (myid == i)
            MPI_Send((renumbering.size() > 0) ? (&renumbering[0]) : nullptr,
                     renumbering.size(),
                     Utilities::MPI::mpi_type_id_for_type<
                       decltype(complete_renumbering[0])>,
                     0,
                     i,
                     MPI_COMM_WORLD);
          else if (myid == 0)
            MPI_Recv((dofs_per_proc[i].n_elements() > 0) ?
                       (&complete_renumbering[offset]) :
                       nullptr,
                     dofs_per_proc[i].n_elements(),
                     Utilities::MPI::mpi_type_id_for_type<
                       decltype(complete_renumbering[0])>,
                     i,
                     i,
                     MPI_COMM_WORLD,
                     MPI_STATUSES_IGNORE);
          offset += dofs_per_proc[i].n_elements();
        }

      if (myid == 0)
        {
          AssertDimension(offset, complete_renumbering.size());
          if (level < tr.n_levels())
            {
              deallog << "On level: " << level << std::endl;
            }
          else
            {
              deallog << "On active cells:" << std::endl;
            }
          for (unsigned int i = 0; i < complete_renumbering.size(); ++i)
            deallog << complete_renumbering[i] << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
