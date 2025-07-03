// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test to check if SolutionTransfer works in parallel with hp::DoFHandler.
// This tests is based on mpi/feindices_transfer.cc


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // ------ setup ------
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(1);

  DoFHandler<dim>       dh(tria);
  hp::FECollection<dim> fe_collection;

  // prepare FECollection with arbitrary number of entries
  const unsigned int max_degree = 1 + Utilities::pow(2, dim);
  for (unsigned int i = 0; i < max_degree; ++i)
    fe_collection.push_back(FE_Q<dim>(max_degree - i));

  typename DoFHandler<dim, dim>::active_cell_iterator cell;
  unsigned int                                        i = 0;

  for (cell = dh.begin_active(); cell != dh.end(); ++cell)
    {
      if (cell->is_locally_owned())
        {
          // set active FE index
          if (!(cell->is_artificial()))
            {
              if (i >= fe_collection.size())
                i = 0;
              cell->set_active_fe_index(i++);
            }

          // set refinement/coarsening flags
          if (cell->id().to_string() == "0_1:0")
            cell->set_refine_flag();
          else if (cell->parent()->id().to_string() ==
                   ((dim == 2) ? "3_0:" : "7_0:"))
            cell->set_coarsen_flag();
        }
    }


  // ----- prepare solution -----
  dh.distribute_dofs(fe_collection);
  IndexSet locally_owned_dofs    = dh.locally_owned_dofs();
  IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dh);

  TrilinosWrappers::MPI::Vector solution;
  solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  for (unsigned int i = 0; i < solution.size(); ++i)
    if (locally_owned_dofs.is_element(i))
      solution(i) = i;

  TrilinosWrappers::MPI::Vector old_solution;
  old_solution.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      MPI_COMM_WORLD);
  old_solution = solution;


  // ----- transfer -----
  SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(dh);

  soltrans.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();

  dh.distribute_dofs(fe_collection);
  locally_owned_dofs    = dh.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dh);

  TrilinosWrappers::MPI::Vector new_solution;
  new_solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  soltrans.interpolate(new_solution);

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
