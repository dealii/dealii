// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test distributed SolutionTransfer with hp-refinement and manually set flags.


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  DoFHandler<dim> dofh(tria);
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(1);

  dofh.distribute_dofs(fes);
  IndexSet locally_owned_dofs = dofh.locally_owned_dofs();
  IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dofh);

  // set up solution
  LinearAlgebra::distributed::Vector<double> solution;
  solution.reinit(locally_owned_dofs,
                  locally_relevant_dofs,
                  dofh.get_mpi_communicator());

  for (unsigned int i = 0; i < solution.size(); ++i)
    if (locally_owned_dofs.is_element(i))
      solution(i) = i;
  solution.update_ghost_values();

  double l1_norm = solution.l1_norm();
  if (Utilities::MPI::this_mpi_process(dofh.get_mpi_communicator()) == 0)
    deallog << "pre  refinement l1=" << l1_norm << std::endl;

  // set refine/coarsen flags manually
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const std::string parent_id_string = cell->parent()->id().to_string();

        if (parent_id_string == "0_0:")
          {
            // parent 0: h refinement
            cell->set_refine_flag();
          }
        else if (parent_id_string == "1_0:")
          {
            // parent 1: h coarsening
            cell->set_coarsen_flag();
          }
        else if (parent_id_string == "2_0:")
          {
            // parent 2: p refinement
            const auto super_fe_index =
              fes.next_in_hierarchy(cell->active_fe_index());
            cell->set_future_fe_index(super_fe_index);
          }
        else if (parent_id_string == "3_0:")
          {
            // parent 3: p coarsening
            const auto sub_fe_index =
              fes.previous_in_hierarchy(cell->active_fe_index());
            cell->set_future_fe_index(sub_fe_index);
          }
        else
          {
            Assert(false, ExcInternalError());
          }
      }

  // initiate refinement and transfer
  SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>> soltrans(
    dofh);
  soltrans.prepare_for_coarsening_and_refinement(solution);

  tria.execute_coarsening_and_refinement();

  dofh.distribute_dofs(fes);
  locally_owned_dofs    = dofh.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dofh);

  solution.reinit(locally_owned_dofs,
                  locally_relevant_dofs,
                  dofh.get_mpi_communicator());
  soltrans.interpolate(solution);

  l1_norm = solution.l1_norm();
  if (Utilities::MPI::this_mpi_process(dofh.get_mpi_communicator()) == 0)
    deallog << "post refinement l1=" << l1_norm << std::endl;

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
