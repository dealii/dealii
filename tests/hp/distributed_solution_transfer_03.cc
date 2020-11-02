// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Test parallel::distributed::SolutionTransfer for FE_Nothing with MPI

#include <deal.II/base/function.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
void
transfer(const MPI_Comm &comm, int test)
{
  parallel::distributed::Triangulation<dim> tria(comm);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  fe.push_back(FE_Nothing<dim>());

  // create a DoFHandler on which we have both cells with FE_Q as well as
  // FE_Nothing
  DoFHandler<dim> dof_handler(tria, true);
  int             i = 0;
  if (test == 0)
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            if (cell->center()[1] < 0.5)
              cell->set_active_fe_index(1);
        }
    }
  else
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              if (i < 2)
                cell->set_active_fe_index(1);
              ++i;
            }
        }
    }
  AffineConstraints<double> cm;
  cm.close();

  dof_handler.distribute_dofs(fe);
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  dof_handler.locally_owned_dofs().print(std::cout);
  locally_relevant_dofs.print(std::cout);

  LinearAlgebra::distributed::Vector<double> solution(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, comm);

  for (unsigned int i = 0; i < solution.local_size(); ++i)
    solution.local_element(i) = i;
  solution.compress(VectorOperation::insert);
  solution.update_ghost_values();

  solution.print(std::cout);

  parallel::distributed::
    SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
      soltrans(dof_handler);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_future_fe_index(0);

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);

  IndexSet new_locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler,
                                          new_locally_relevant_dofs);
  LinearAlgebra::distributed::Vector<double> interpolated_solution(
    dof_handler.locally_owned_dofs(), new_locally_relevant_dofs, comm);

  soltrans.interpolate(interpolated_solution);

  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "-------------------------" << std::endl;
  solution.print(std::cout);

  // if (test == 0)
  //   {
  //     std::vector<double> ref = {0., 1., 2., 3., 4., 5., 0., 0., 0.};
  //     for (unsigned int i = 0; i < ref.size(); ++i)
  //       {
  //         AssertThrow(std::abs(solution[i] - ref[i]) < 1e-12,
  //                     ExcInternalError());
  //       }
  //   }
  // else
  //   {
  //     std::vector<double> ref = {0., 0., 0., 1., 0., 4., 2., 3., 5.};
  //     for (unsigned int i = 0; i < ref.size(); ++i)
  //       AssertThrow(std::abs(solution[i] - ref[i]) < 1e-12,
  //       ExcInternalError());
  //   }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const MPI_Comm comm = MPI_COMM_WORLD;
  initlog();

  deallog.push("2D solution transfer");
  transfer<2>(comm, 0);
  //  transfer<2>(comm, 1);
  deallog << "OK" << std::endl;
  deallog.pop();
}
