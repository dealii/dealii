// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// save and load a triangulation with one solution vector

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
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  std::string filename = "dat";
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

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

    FE_Q<dim>       fe(1);
    DoFHandler<dim> dh(tr);

    dh.distribute_dofs(fe);

    const IndexSet &locally_owned_dofs = dh.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dh);

    PETScWrappers::MPI::Vector x(locally_owned_dofs, MPI_COMM_WORLD);
    PETScWrappers::MPI::Vector solution(locally_owned_dofs,
                                        locally_relevant_dofs,
                                        MPI_COMM_WORLD);

    SolutionTransfer<dim, PETScWrappers::MPI::Vector> soltrans(dh);

    for (unsigned int i = 0; i < locally_owned_dofs.n_elements(); ++i)
      {
        unsigned int idx = locally_owned_dofs.nth_index_in_set(i);
        x(idx)           = idx;

        //  std::cout << '[' << idx << ']' << ' ' << solution(idx) << std::endl;
      }

    x.compress(VectorOperation::insert);
    solution = x;

    soltrans.prepare_for_serialization(solution);

    tr.save(filename.c_str());

    if (myid == 0)
      {
        deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
        deallog << "cells(0) = " << tr.n_active_cells() << std::endl;
      }
    deallog << "Checksum: " << tr.get_checksum() << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.load(filename.c_str());
    FE_Q<dim>       fe(1);
    DoFHandler<dim> dh(tr);

    dh.distribute_dofs(fe);

    const IndexSet &locally_owned_dofs = dh.locally_owned_dofs();
    const IndexSet  locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dh);

    PETScWrappers::MPI::Vector solution(locally_owned_dofs, MPI_COMM_WORLD);
    SolutionTransfer<dim, PETScWrappers::MPI::Vector> soltrans(dh);
    solution = 2;
    soltrans.deserialize(solution);

    for (unsigned int i = 0; i < locally_owned_dofs.n_elements(); ++i)
      {
        unsigned int idx = locally_owned_dofs.nth_index_in_set(i);
        // std::cout << '[' << idx << ']' << ' ' << solution(idx) << std::endl;
        Assert(idx == get_real_assert_zero_imag(solution(idx)),
               ExcInternalError());
      }

    double norm = solution.l1_norm();

    if (myid == 0)
      {
        deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
        deallog << "cells(0) = " << tr.n_active_cells() << std::endl;
      }
    deallog << "Checksum: " << tr.get_checksum() << std::endl;
    deallog << "sum: " << norm << std::endl;
  }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
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
    }
  else
    test<2>();
}
