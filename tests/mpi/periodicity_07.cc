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

// a variation of periodicity_06 that used to trigger
// AffineConstraints::is_consistent_in_parallel() on 13 mpi tasks.

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
test(const unsigned numRefinementLevels = 2)
{
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);

  ConditionalOStream pcout(std::cout, this_mpi_process == 0);

  const double                                      L = 20;
  dealii::parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator);
  GridGenerator::hyper_cube(triangulation, -L, L, /*colorize*/ false);

  // mark faces
  for (auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      {
        const Point<dim> &face_center = cell->face(f)->center();
        if (cell->face(f)->at_boundary())
          {
            unsigned int counter = 1;
            for (unsigned int d = 0; d < dim; ++d)
              {
                if (std::abs(face_center[d] - L) < 1.0e-5)
                  cell->face(f)->set_boundary_id(counter);
                ++counter;
                if (std::abs(face_center[d] + L) < 1.0e-5)
                  cell->face(f)->set_boundary_id(counter);
                ++counter;
              }
          }
      }

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  for (int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(triangulation,
                                      /*b_id1*/ 2 * d + 1,
                                      /*b_id2*/ 2 * d + 2,
                                      /*direction*/ d,
                                      periodicity_vector);

  triangulation.add_periodicity(periodicity_vector);

  // refine mesh
  triangulation.refine_global(1);

  Point<dim> corner;
  for (unsigned int d = 0; d < dim; ++d)
    corner[d] = -L;

  MappingQ1<dim> mapping;
  for (unsigned int ilevel = 0; ilevel < numRefinementLevels; ++ilevel)
    {
      // pick an corner cell and refine
      for (auto &cell : triangulation.active_cell_iterators())
        {
          try
            {
              const Point<dim> p_cell =
                mapping.transform_real_to_unit_cell(cell, corner);
              const double dist =
                GeometryInfo<dim>::distance_to_unit_cell(p_cell);

              if (dist < 1e-08)
                cell->set_refine_flag();
            }
          catch (const typename MappingQ1<dim>::ExcTransformationFailed &)
            {}
        }
      triangulation.execute_coarsening_and_refinement();
    }

  pcout << "number of elements: " << triangulation.n_global_active_cells()
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

  std::vector<
    GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
    periodicity_vectorDof;
  for (int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(dof_handler,
                                      /*b_id1*/ 2 * d + 1,
                                      /*b_id2*/ 2 * d + 2,
                                      /*direction*/ d,
                                      periodicity_vectorDof);

  DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorDof,
                                                   constraints);

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
  test<3>(4);
}
