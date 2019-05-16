/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * A 2d example of smoothing of a distorted grid using the Mesquite interface.
 * The grid is locally refined, and so it has hanging nodes.
 * This test is similar to mesh_interface_tria_07, except that we now test
 * a parallel triangulation and no refinement cycles are considered.
 *
 * - The vertices of the triangulation are moved to a distorted position
 * - Mesh adaption is performed
 *     - A predefined smoother is used
 *     - The boundary vertices are considered fixed
 *     - We rely on the automatic enforcement of hanging vertex constraints at
 *       the end of MesquiteMeshInterface::execute()
 * - The vertices of the triangulation are updated to their optimial locations
 */


#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>

#include <fstream>

#include "../tests.h"
#include "mesquite_utilities.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  const MPI_Comm &   mpi_communicator = MPI_COMM_WORLD;
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  ConditionalOStream pcout(std::cout, this_mpi_process == 0);

  const unsigned int dim                 = 2;
  const unsigned int n_elements_per_side = 4;
  const unsigned int n_refinement_cycles = 0;

  // Make a distored grid
  parallel::shared::Triangulation<dim> tria(
    mpi_communicator, Triangulation<dim>::maximum_smoothing);
  // There's a bug in distort random, so for now we use another distorted grid
  const bool use_distort_random = false;
  if (use_distort_random)
    {
      std::vector<unsigned int> repetitions(dim, n_elements_per_side);
      GridGenerator::subdivided_hyper_rectangle(tria,
                                                repetitions,
                                                Point<2>(0.0, 0.0),
                                                Point<2>(1.0, 1.0));
      GridTools::distort_random(0.3, tria, true);
    }
  else
    {
      // This is the same grid as introduced in mesh_interface_tria_04
      std::vector<unsigned int> repetitions(dim);
      repetitions[0] = 14;
      repetitions[1] = 4;
      GridGenerator::subdivided_hyper_rectangle(tria,
                                                repetitions,
                                                Point<2>(0.0, 0.0),
                                                Point<2>(14.0, 2.0));
      //    GridTools::transform(TransformationSinusoid<dim>(false), tria);

      // Prescribe the mesh motion
      const unsigned int   first_u_dof = 0;
      const Vector<double> eulerian_vertex_map_exact =
        GridTools::get_eulerian_vertex_positions(
          tria, TransformationSinusoid<dim>(false), first_u_dof);
      const Vector<double> eulerian_vertex_map_ptrb =
        GridTools::get_eulerian_vertex_positions(
          tria, TransformationSinusoid<dim>(true), first_u_dof);

      // Undo the perturbations of the boundary vertices
      Vector<double> eulerian_vertex_map(dim * tria.n_vertices());
      fix_eulerian_boundary_vertex_map(eulerian_vertex_map,
                                       tria,
                                       eulerian_vertex_map_exact,
                                       eulerian_vertex_map_ptrb);
      GridTools::move_triangulation_vertices(tria, eulerian_vertex_map);
    }

  for (unsigned int c = 0; c < n_refinement_cycles; ++c)
    {
      tria.begin_active(c)->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  // Output the distorted grid
  pcout << "Distorted mesh" << std::endl;
  output_mesh(tria, "grid.01." + Utilities::int_to_string(this_mpi_process));
  //  output_mesh_vertex_numbering (tria, "grid.01." +
  //  Utilities::int_to_string(this_mpi_process));

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    tria, true /*fix_all_boundary_vertices*/);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::laplace);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  pcout << "Smoothed mesh" << std::endl;
  output_mesh(tria,
              "grid.02." + Utilities::int_to_string(this_mpi_process),
              true /*write_to_deallog*/);

  pcout << "OK" << std::endl;
}
