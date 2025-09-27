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



// check that inserting particles works:
// - in parallel when one processor is not inserting any particles
// - and when all cells are locally owned on one processor

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(const unsigned int refinement)
{
  {
    parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

    GridGenerator::hyper_cube(tr);
    tr.refine_global(refinement);
    MappingQ<dim, spacedim> mapping(1);

    deallog << "Refinement: " << refinement << std::endl;

    deallog << "Global active cells: " << tr.n_global_active_cells()
            << std::endl;
    deallog << "Locally owned active cells: "
            << tr.n_locally_owned_active_cells() << std::endl;

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<dim>> positions(1, Point<dim>());
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
      positions.clear();

    particle_handler.insert_particles(positions);

    deallog << "Global particles: " << particle_handler.n_global_particles()
            << std::endl;
    deallog << "Locally owned particles: "
            << particle_handler.n_locally_owned_particles() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  deallog.push("2d/2d");
  test<2, 2>(1);
  deallog.pop();
  deallog << "---" << std::endl;
  deallog.push("2d/2d");
  test<2, 2>(2);
  deallog.pop();
}
