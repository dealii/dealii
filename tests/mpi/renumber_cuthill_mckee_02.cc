// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that DofRenumbering::Cuthill_McKee works in parallel also when
// a set of starting indices is given.


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

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

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
  std::vector<types::global_dof_index> renumbering(
    dofh.locally_owned_dofs().n_elements());
  std::vector<types::global_dof_index> starting_indices;
  starting_indices.push_back(dofh.locally_owned_dofs().nth_index_in_set(0));
  DoFRenumbering::compute_Cuthill_McKee(
    renumbering, dofh, false, false, starting_indices);

  // send everything to processor 0 for output
  std::vector<types::global_dof_index> complete_renumbering(dofh.n_dofs());
  std::copy(renumbering.begin(),
            renumbering.end(),
            complete_renumbering.begin());
  unsigned int                offset = renumbering.size();
  const std::vector<IndexSet> locally_owned_dofs_per_processor =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, dofh.locally_owned_dofs());
  for (unsigned int i = 1; i < nprocs; ++i)
    {
      if (myid == i)
        MPI_Send(&renumbering[0],
                 renumbering.size(),
                 Utilities::MPI::mpi_type_id_for_type<
                   decltype(complete_renumbering[0])>,
                 0,
                 i,
                 MPI_COMM_WORLD);
      else if (myid == 0)
        MPI_Recv(&complete_renumbering[offset],
                 locally_owned_dofs_per_processor[i].n_elements(),
                 Utilities::MPI::mpi_type_id_for_type<
                   decltype(complete_renumbering[0])>,
                 i,
                 i,
                 MPI_COMM_WORLD,
                 MPI_STATUSES_IGNORE);
      offset += locally_owned_dofs_per_processor[i].n_elements();
    }

  if (myid == 0)
    {
      AssertDimension(offset, complete_renumbering.size());
      for (unsigned int i = 0; i < complete_renumbering.size(); ++i)
        deallog << complete_renumbering[i] << std::endl;
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
