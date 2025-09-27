// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test if sort_particles_into_subdomains_and_cells() works for a particle
// located at a vertex.

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::subdivided_hyper_cube(tr, 20, 0, 4);

  MappingQ<dim, spacedim> mapping(3);

  // both processes create a particle handler with a particle in a
  // position that is in a ghost cell to the other process
  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  particle_handler.signals.particle_lost.connect(
    [](const typename Particles::ParticleIterator<dim>         &particle,
       const typename Triangulation<dim>::active_cell_iterator &cell) {
      AssertThrow(false, ExcMessage("Particle lost."));
    });

  std::vector<Point<spacedim>> position(1);

  if (Utilities::MPI::this_mpi_process(tr.get_mpi_communicator()) == 0)
    for (unsigned int i = 0; i < dim; ++i)
      position[0][i] = 2;

  auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    tr, IteratorFilters::LocallyOwnedCell());

  auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  particle_handler.insert_global_particles(position, global_bounding_boxes);

  // Reverse the refinement and check again
  for (auto cell = tr.begin_active(); cell != tr.end(); ++cell)
    {
      if (cell->is_locally_owned())
        if (cell->center().distance(position[0]) <= 1)
          cell->set_refine_flag();
    }

  particle_handler.prepare_for_coarsening_and_refinement();
  tr.prepare_coarsening_and_refinement();
  tr.execute_coarsening_and_refinement();
  particle_handler.unpack_after_coarsening_and_refinement();
  particle_handler.sort_particles_into_subdomains_and_cells();

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  test<2, 2>();
}
