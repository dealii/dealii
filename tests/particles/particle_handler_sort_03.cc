// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test that ParticleHandler::sort_particles_into_subdomains_and_cells() does
// not silently drop a particle that moves into a cell that is artificial on
// the process that owns the particle. Such a particle used to be scheduled to
// be sent to numbers::artificial_subdomain_id, which is not a real rank: in
// debug mode this triggered an assertion, in release mode the particle
// disappeared without the particle_lost signal ever being emitted. Artificial
// cells are no valid target, so the particle has to be reported as lost.
//
// The setup is two coarse cells merged into [-1,1]x[0,1] and refined once,
// i.e. eight cells. On two processes rank 0 owns the four cells with x < 0
// and its ghost layer only covers 0 < x < 0.5, so that the two cells with
// x > 0.5 are artificial on rank 0. A single particle owned by rank 0 is
// moved from (-0.75, 0.25) into one of those artificial cells in one step.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

void
test()
{
  const unsigned int dim = 2;

  const MPI_Comm     comm    = MPI_COMM_WORLD;
  const unsigned int my_rank = Utilities::MPI::this_mpi_process(comm);

  const Point<dim> start_location(-0.75, 0.25);
  const Point<dim> end_location(0.75, 0.25);

  parallel::distributed::Triangulation<dim> tria(comm);
  GridGenerator::subdivided_hyper_rectangle(
    tria, {2, 1}, Point<dim>(-1., 0.), Point<dim>(1., 1.), true);
  tria.refine_global(1);

  const MappingQ1<dim>            mapping;
  Particles::ParticleHandler<dim> particle_handler(tria, mapping);

  // Count how often the library tells us that a particle got lost. Before the
  // fix this counter stayed at zero even though the particle disappeared.
  unsigned int n_lost_particles = 0;
  particle_handler.signals.particle_lost.connect(
    [&n_lost_particles](
      const typename Particles::ParticleIterator<dim>         &particle,
      const typename Triangulation<dim>::active_cell_iterator &cell) {
      ++n_lost_particles;
      deallog << "Particle <" << particle->get_id() << "> lost in cell "
              << cell->id() << std::endl;
    });

  // Insert the one and only particle, owned by rank 0.
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    tria, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(comm, my_bounding_box);

  std::vector<Point<dim>> positions;
  if (my_rank == 0)
    positions.push_back(start_location);
  particle_handler.insert_global_particles(positions, global_bounding_boxes);

  deallog << "Particles before sort: " << particle_handler.n_global_particles()
          << std::endl;

  // The single displacement of the single particle. Its new location lies in
  // a cell that is artificial on rank 0.
  for (auto &particle : particle_handler)
    particle.get_location() = end_location;

  particle_handler.sort_particles_into_subdomains_and_cells();

  deallog << "Particles after sort:  " << particle_handler.n_global_particles()
          << std::endl;
  deallog << "Lost particles:        "
          << Utilities::MPI::sum(n_lost_particles, comm) << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  test();
}
