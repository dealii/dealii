// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test adaptive refinement with periodic boundary conditions and
// mesh_reconstruction_after_repartitioning

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::none,
    typename parallel::distributed::Triangulation<dim>::Settings(
      parallel::distributed::Triangulation<
        dim>::mesh_reconstruction_after_repartitioning |
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy |
      parallel::distributed::Triangulation<
        dim>::communicate_vertices_to_p4est));

  GridGenerator::subdivided_hyper_cube(tr, 4, -1.0, 3.0, true);

  // make x periodic
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tr, 0, 1, 0, periodic_faces);

  tr.add_periodicity(periodic_faces);

  for (unsigned int i = 0; i < 6; ++i)
    {
      deallog << "Checksum: " << tr.get_checksum() << '\n'
              << "n_active_cells: " << tr.n_active_cells() << '\n'
              << "n_levels: " << tr.n_levels() << std::endl;

      Point<dim> p(-0.25, 0);

      for (auto &cell : tr.active_cell_iterators())
        {
          if (p.distance(cell->center()) < 0.75)
            cell->set_refine_flag();
        }

      tr.execute_coarsening_and_refinement();
    }
  tr.write_mesh_vtk("mesh");
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  test<2>();
}
