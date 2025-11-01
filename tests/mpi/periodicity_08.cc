/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2019 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Document a hang in Triangulation::prepare_coarsening_and_refinement() when
// periodic boundaries and mesh smoothing is used.

// The (now fixed bug) was caused by eliminate_refined_boundary_islands also
// incorrectly acting on periodic boundaries. This test originally comes from
// ASPECT.

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

  const double L = 20;

  dealii::parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::limit_level_difference_at_vertices |
      Triangulation<dim>::eliminate_unrefined_islands |
      Triangulation<dim>::eliminate_refined_inner_islands |
      Triangulation<dim>::eliminate_refined_boundary_islands |
      Triangulation<dim>::do_not_produce_unrefined_islands));
  GridGenerator::hyper_cube(triangulation, 0, 1, /*colorize*/ true);



  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  if (true)
    for (int d = 0; d < dim; ++d)
      GridTools::collect_periodic_faces(triangulation,
                                        /*b_id1*/ 2 * d,
                                        /*b_id2*/ 2 * d + 1,
                                        /*direction*/ d,
                                        periodicity_vector);

  triangulation.add_periodicity(periodicity_vector);

  // refine mesh
  triangulation.refine_global(3);

  // mark all cells except these for coarsening:
  std::vector<Point<dim>> points = {{312500, 93750},
                                    {312500, 156250},
                                    {62500, 343750},
                                    {312500, 281250},
                                    {312500, 343750},
                                    {62500, 406250},
                                    {62500, 468750}};

  // rescale to [0,1]^2 for this to work:
  for (auto &p : points)
    {
      p[0] /= 1e6;
      p[1] /= 5e5;
    }

  for (auto &cell : triangulation.active_cell_iterators())
    {
      cell->set_coarsen_flag();
      for (auto &p : points)
        if (cell->point_inside(p))
          {
            cell->clear_coarsen_flag();
            cell->set_refine_flag();
          }
    }
  triangulation.execute_coarsening_and_refinement();

  deallog << "number of elements: " << triangulation.n_global_active_cells()
          << std::endl;

  // create dof_handler
  FESystem<dim>   FE(FE_Q<dim>(QGaussLobatto<1>(2)), 1);
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

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  const IndexSet locally_active_dofs =
    DoFTools::extract_locally_active_dofs(dof_handler);

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
  constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);

  for (int d = 0; d < dim; ++d)
    DoFTools::make_periodicity_constraints(
      dof_handler, 2 * d, 2 * d + 1, d, constraints);

  const bool consistent =
    constraints.is_consistent_in_parallel(locally_owned_dofs,
                                          locally_active_dofs,
                                          mpi_communicator,
                                          /*verbose*/ true);

  deallog << "Periodicity constraints are consistent in parallel: "
          << consistent << std::endl;

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  const bool hanging_consistent =
    constraints.is_consistent_in_parallel(locally_owned_dofs,
                                          locally_active_dofs,
                                          mpi_communicator);

  deallog << "Hanging nodes constraints are consistent in parallel: "
          << hanging_consistent << std::endl;

  constraints.close();
  deallog << "OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    log;
  test<2>();
}
