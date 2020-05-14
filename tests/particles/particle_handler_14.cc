// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// check variable size data transfer for parallel distributed applications
// if one cell does not contain any particle at all


#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  MappingQ<dim, spacedim> mapping(1);

  // setup triangulation
  const unsigned int n_cells = 2;

  std::vector<unsigned int> rep(dim, 1);
  rep[0] = n_cells;
  Point<dim> p1, p2;
  for (unsigned int d = 0; d < dim; ++d)
    {
      p1[d] = 0;
      p2[d] = (d == 0) ? n_cells : 1;
    }
  GridGenerator::subdivided_hyper_rectangle(tr, rep, p1, p2);

  // setup one particle in first cell
  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  Particles::Particle<dim, spacedim> particle(Point<spacedim>(),
                                              Point<dim>(),
                                              0);

  particle_handler.insert_particle(particle, tr.begin_active());
  particle_handler.update_cached_numbers();

  // initiate data transfer
  tr.signals.pre_distributed_repartition.connect([&particle_handler]() {
    particle_handler.register_store_callback_function();
  });

  tr.signals.post_distributed_repartition.connect([&particle_handler]() {
    particle_handler.register_load_callback_function(false);
  });

  tr.repartition();

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    initlog();

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
