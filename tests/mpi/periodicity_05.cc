/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2017 - 2020 by the deal.II authors
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


// This testcase previously failed due to missing to update the NumberCache in
// add_periodicity(). This resulted in ghost_owners to be wrong and not sending
// all required messages in GridTools::exchange_cell_data_to_ghosts.

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  mpi_initlog();
  const double                            L = 7.6 / 2.0;
  parallel::distributed::Triangulation<3> triangulation(MPI_COMM_WORLD);
  GridIn<3>                               gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f(SOURCE_DIR "/periodicity_05.inp");
  gridin.read_ucd(f);

  typename parallel::distributed::Triangulation<3>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (; cell != endc; ++cell)
    {
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          const Point<3> face_center = cell->face(f)->center();
          if (cell->face(f)->at_boundary())
            {
              if (std::abs(face_center[0] + (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(1);
              else if (std::abs(face_center[0] - (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(2);
              else if (std::abs(face_center[1] + (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(3);
              else if (std::abs(face_center[1] - (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(4);
              else if (std::abs(face_center[2] + (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(5);
              else if (std::abs(face_center[2] - (L)) < 1.0e-4)
                cell->face(f)->set_boundary_id(6);
            }
        }
    }

  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<3>::cell_iterator>>
    periodicity_vector;
  for (int i = 0; i < 3; ++i)
    GridTools::collect_periodic_faces(triangulation,
                                      /*b_id1*/ 2 * i + 1,
                                      /*b_id2*/ 2 * i + 2,
                                      /*direction*/ i,
                                      periodicity_vector);
  triangulation.add_periodicity(periodicity_vector);

  FE_Q<3>       fe(1);
  DoFHandler<3> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  deallog << "OK" << std::endl;
}
