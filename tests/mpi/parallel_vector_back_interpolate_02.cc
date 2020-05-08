// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


// check FETools::back_interpolate on Trilinos vectors

#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int dim = 2;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FESystem<dim>   fe1(FE_Q<dim>(1), 1, FE_Q<dim>(1), 1);
  FESystem<dim>   fe2(FE_Q<dim>(2), 1, FE_Q<dim>(2), 1);
  DoFHandler<dim> dof1(tria), dof2(tria);
  dof1.distribute_dofs(fe1);
  dof2.distribute_dofs(fe2);

  DoFRenumbering::component_wise(dof1);
  DoFRenumbering::component_wise(dof2);

  IndexSet locally_owned_dofs1 = dof1.locally_owned_dofs();
  IndexSet locally_relevant_dofs1;
  DoFTools::extract_locally_relevant_dofs(dof1, locally_relevant_dofs1);

  std::vector<IndexSet> owned_partitioning1;
  std::vector<IndexSet> relevant_partitioning1;
  {
    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof1);
    const unsigned int n1 = dofs_per_block[0], n2 = dofs_per_block[1];

    owned_partitioning1.push_back(locally_owned_dofs1.get_view(0, n1));
    owned_partitioning1.push_back(locally_owned_dofs1.get_view(n1, n1 + n2));

    relevant_partitioning1.push_back(locally_relevant_dofs1.get_view(0, n1));
    relevant_partitioning1.push_back(
      locally_relevant_dofs1.get_view(n1, n1 + n2));
  }

  AffineConstraints<double> c1, c2;
  DoFTools::make_hanging_node_constraints(dof1, c1);
  c1.close();
  DoFTools::make_hanging_node_constraints(dof2, c2);
  c2.close();

  {
    TrilinosWrappers::MPI::Vector v1_distributed(locally_owned_dofs1,
                                                 MPI_COMM_WORLD);
    TrilinosWrappers::MPI::Vector v1(locally_owned_dofs1,
                                     locally_relevant_dofs1,
                                     MPI_COMM_WORLD);
    TrilinosWrappers::MPI::Vector v1_interpolated(v1_distributed);

    for (const auto &el : locally_owned_dofs1)
      v1_distributed(el) = random_value();
    c1.distribute(v1_distributed);
    v1 = v1_distributed;

    FETools::back_interpolate(dof1, c1, v1, dof2, c2, v1_interpolated);
    for (const auto &el : locally_owned_dofs1)
      {
        if (std::abs(v1_interpolated(el) - v1(el)) > 1.e-10)
          {
            std::cout << v1_interpolated(el) << " should be " << v1(el)
                      << std::endl;
            AssertThrow(false, ExcInternalError());
          }
      }
    deallog << "TrilinosWrappers::MPI::Vector: OK" << std::endl;
  }
  {
    TrilinosWrappers::MPI::BlockVector v1_distributed(owned_partitioning1,
                                                      MPI_COMM_WORLD);
    TrilinosWrappers::MPI::BlockVector v1(owned_partitioning1,
                                          relevant_partitioning1,
                                          MPI_COMM_WORLD);
    TrilinosWrappers::MPI::BlockVector v1_interpolated(v1_distributed);

    for (const auto &el : locally_owned_dofs1)
      v1_distributed(el) = random_value();
    c1.distribute(v1_distributed);
    v1 = v1_distributed;

    FETools::back_interpolate(dof1, c1, v1, dof2, c2, v1_interpolated);
    for (const auto &el : locally_owned_dofs1)
      {
        if (std::abs(v1_interpolated(el) - v1(el)) > 1.e-10)
          {
            std::cout << v1_interpolated(el) << " should be " << v1(el)
                      << std::endl;
            AssertThrow(false, ExcInternalError());
          }
      }
    deallog << "TrilinosWrappers::MPI::BlockVector: OK" << std::endl;
  }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
