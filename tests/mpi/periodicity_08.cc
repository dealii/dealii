/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// A test that triggers invalid contraints with 10 MPI ranks

/*

-------------------------------------------------------
An error occurred in line <408> of file
<../include/deal.II/lac/affine_constraints.templates.h> in function void
dealii::AffineConstraints<number>::close() [with number = double] The violated
condition was: dof_index != line.index Additional information: Cycle in
constraints detected!
--------------------------------------------------------


 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

template <int dim>
void
test()
{
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);

  ConditionalOStream pcout(std::cout, this_mpi_process == 0);

  dealii::parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator);
  std::vector<unsigned int> repetitions = {9, 3, 4};
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            Point<dim>(0., 0., 0.),
                                            Point<dim>(9.0, 3.0, 4.0),
                                            /*colorize*/ true);

  // make y periodic:
  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  GridTools::collect_periodic_faces(triangulation,
                                    /*b_id1*/ 2 * 1,
                                    /*b_id2*/ 2 * 1 + 1,
                                    /*direction*/ 1,
                                    periodicity_vector);

  triangulation.add_periodicity(periodicity_vector);


  Point<dim> p(8.5, 0.5, 0.5);

  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->center().distance(p) > 0.1)
        cell->set_refine_flag();
    }

  triangulation.execute_coarsening_and_refinement();

  pcout << "number of elements: " << triangulation.n_global_active_cells()
        << std::endl;

  // create dof_handler
  FESystem<dim>   FE(FE_Q<dim>(2), 3, FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE);

  // write mesh for visualization
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
  data_out.build_patches();
  data_out.write_vtu_in_parallel(std::string("mesh.vtu").c_str(),
                                 mpi_communicator);

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  IndexSet locally_active_dofs;
  DoFTools::extract_locally_active_dofs(dof_handler, locally_active_dofs);

  const std::vector<IndexSet> locally_owned_dofs =
    Utilities::MPI::all_gather(MPI_COMM_WORLD,
                               dof_handler.locally_owned_dofs());

  std::map<types::global_dof_index, Point<dim>> supportPoints;
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),
                                       dof_handler,
                                       supportPoints);

  /// creating combined hanging node and periodic constraint matrix
  AffineConstraints<double> constraints;
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);

  DoFTools::make_periodicity_constraints(
    dof_handler, 2 * 1, 2 * 1 + 1, 1, constraints);

  const bool consistent =
    constraints.is_consistent_in_parallel(locally_owned_dofs,
                                          locally_active_dofs,
                                          mpi_communicator,
                                          /*verbose*/ true);

  pcout << "Periodicity constraints are consistent in parallel: " << consistent
        << std::endl;

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  const bool hanging_consistent =
    constraints.is_consistent_in_parallel(locally_owned_dofs,
                                          locally_active_dofs,
                                          mpi_communicator);

  pcout << "Hanging nodes constraints are consistent in parallel: "
        << hanging_consistent << std::endl;

  constraints.close();
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  test<3>();
}
