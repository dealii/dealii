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

// a unit test from https://github.com/dealii/dealii/issues/7053
// we used to build incompatible constraints:
// === Process 0
// Constraints:
// ...
//     16 32:  1
// ...
//     52 32:  1
// ...
// Owned DoFs:
// {[0,93]}
// Ghost DoFs:
// {[0,143]}
// Coordinates:
// 52@-20 20 -10
// 32@20 -20 -10
//
// vs
//
// === Process 1
// Constraints:
// ...
//     32 16:  1
// ...
//     52 16:  1
// ...
// Owned DoFs:
// {[94,143]}
// Ghost DoFs:
// {[0,88], [94,143]}
// Coordinates:
// 52@-20 20 -10
// 16@-20 -20 -10

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
  for (unsigned int ilevel = 0; ilevel < numRefinementLevels; ilevel++)
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

  if (this_mpi_process == 0)
    deallog << "number of elements: " << triangulation.n_global_active_cells()
            << '\n';

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

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  std::map<types::global_dof_index, Point<dim>> supportPoints;
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),
                                       dof_handler,
                                       supportPoints);

  /// creating combined hanging node and periodic constraint matrix
  AffineConstraints<double> constraints;
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  std::vector<
    GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
    periodicity_vectorDof;
  for (int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(dof_handler,
                                      /*b_id1*/ 2 * d + 1,
                                      /*b_id2*/ 2 * d + 2,
                                      /*direction*/ d,
                                      periodicity_vectorDof);

  DoFTools::make_periodicity_constraints<DoFHandler<dim>>(periodicity_vectorDof,
                                                          constraints);
  constraints.close();

  const std::vector<IndexSet> &locally_owned_dofs =
    Utilities::MPI::all_gather(MPI_COMM_WORLD,
                               dof_handler.locally_owned_dofs());
  IndexSet locally_active_dofs;
  DoFTools::extract_locally_active_dofs(dof_handler, locally_active_dofs);
  AssertThrow(constraints.is_consistent_in_parallel(locally_owned_dofs,
                                                    locally_active_dofs,
                                                    mpi_communicator,
                                                    /*verbose*/ true),
              ExcInternalError());

  deallog << "=== Process " << this_mpi_process << std::endl
          << "Constraints:" << std::endl;
  constraints.print(deallog.get_file_stream());
  deallog << "Owned DoFs:" << std::endl;
  dof_handler.locally_owned_dofs().print(deallog);
  deallog << "Ghost DoFs:" << std::endl;
  locally_relevant_dofs.print(deallog);
  // we have issues with constraints for u_52 = ....
  const unsigned int i = 52;
  if (locally_relevant_dofs.is_element(i) && constraints.is_constrained(i))
    {
      deallog << "Coordinates:" << std::endl;
      deallog << i << "@" << supportPoints[i] << std::endl;
      for (auto c : *constraints.get_constraint_entries(i))
        deallog << c.first << "@" << supportPoints[c.first] << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    mpi_log;
  test<3>();
}
