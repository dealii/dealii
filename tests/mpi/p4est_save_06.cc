// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// save and load a triangulation with a different number of cpus
// using variable size data transfer with an associated hp::DoFHandler

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  MPI_Comm     com_all = MPI_COMM_WORLD;
  MPI_Comm     com_small;

  // split the communicator in proc 0,1,2 and 3,4
  MPI_Comm_split(com_all, (myid < 3) ? 0 : 1, myid, &com_small);

  // write with small com
  if (myid < 3)
    {
      deallog << "writing with " << Utilities::MPI::n_mpi_processes(com_small)
              << std::endl;

      parallel::distributed::Triangulation<dim> tr(com_small);
      GridGenerator::hyper_cube(tr);
      tr.refine_global(2);
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

      DoFHandler<dim>       dh(tr);
      hp::FECollection<dim> fe_collection;

      // prepare FECollection with arbitrary number of entries
      const unsigned int max_degree = 1 + Utilities::pow(2, dim);
      for (unsigned int i = 0; i < max_degree; ++i)
        fe_collection.push_back(FE_Q<dim>(max_degree - i));

      dh.distribute_dofs(fe_collection);

      const IndexSet &locally_owned_dofs = dh.locally_owned_dofs();
      const IndexSet  locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dh);

      PETScWrappers::MPI::Vector x(locally_owned_dofs, com_small);
      PETScWrappers::MPI::Vector rel_x(locally_owned_dofs,
                                       locally_relevant_dofs,
                                       com_small);

      SolutionTransfer<dim, PETScWrappers::MPI::Vector> soltrans(dh);

      for (unsigned int i = 0; i < locally_owned_dofs.n_elements(); ++i)
        {
          unsigned int idx = locally_owned_dofs.nth_index_in_set(i);
          x(idx)           = idx;
        }


      x.compress(VectorOperation::insert);
      rel_x = x;

      dh.prepare_for_serialization_of_active_fe_indices();
      soltrans.prepare_for_serialization(rel_x);

      tr.save("file");
      deallog << "#cells: " << tr.n_global_active_cells() << std::endl
              << "norm: " << x.l2_norm() << std::endl
              << "Checksum: " << tr.get_checksum() << std::endl;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "reading with " << Utilities::MPI::n_mpi_processes(com_all)
          << std::endl;

  {
    parallel::distributed::Triangulation<dim> tr(com_all);

    GridGenerator::hyper_cube(tr);
    tr.load("file");

    DoFHandler<dim>       dh(tr);
    hp::FECollection<dim> fe_collection;

    // prepare FECollection with arbitrary number of entries
    const unsigned int max_degree = 1 + Utilities::pow(2, dim);
    for (unsigned int i = 0; i < max_degree; ++i)
      fe_collection.push_back(FE_Q<dim>(max_degree - i));


    dh.deserialize_active_fe_indices();
    dh.distribute_dofs(fe_collection);

    const IndexSet &locally_owned_dofs = dh.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dh);

    PETScWrappers::MPI::Vector solution(locally_owned_dofs, com_all);
    solution = PetscScalar();

    SolutionTransfer<dim, PETScWrappers::MPI::Vector> soltrans(dh);

    soltrans.deserialize(solution);

    deallog << "#cells: " << tr.n_global_active_cells() << std::endl
            << "norm: " << solution.l2_norm() << std::endl
            << "Checksum: " << tr.get_checksum() << std::endl;
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
