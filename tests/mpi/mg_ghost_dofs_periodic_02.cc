// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// test distribute_mg_dofs and the ghost layer for periodic boundary
// conditions
// same as mg_ghst_dofs_periodic_01 but calling add_periodicity twice and only
// in 3D

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Assert(dim == 3, ExcNotImplemented());
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria, 0, 1);
  for (unsigned int face = 2; face < GeometryInfo<dim>::faces_per_cell; ++face)
    tria.begin()->face(face)->set_all_boundary_ids(face);

  for (unsigned int d = 1; d < dim; ++d)
    {
      std::vector<
        GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
        periodic_faces;
      GridTools::collect_periodic_faces(
        tria, 2 * d, 2 * d + 1, d, periodic_faces);
      tria.add_periodicity(periodic_faces);
    }

  tria.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      deallog << "Level " << level << std::endl;
      for (typename DoFHandler<dim>::cell_iterator cell =
             dof_handler.begin(level);
           cell != dof_handler.end(level);
           ++cell)
        if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            deallog << "Cell with center: " << cell->center() << ", owned by "
                    << cell->level_subdomain_id() << ": ";
            cell->get_mg_dof_indices(dof_indices);
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              deallog << dof_indices[i] << " ";
            deallog << std::endl;
          }
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPILogInitAll log;
  deallog << std::setprecision(4);

  try
    {
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
