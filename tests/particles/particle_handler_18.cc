// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like particle_handler_08, but check that the properties of particles are
// correctly copied in the ParticleHandler::copy_from function.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tr;

  GridGenerator::hyper_cube(tr);
  MappingQ<dim, spacedim>                   mapping(1);
  const unsigned int                        n_properties = spacedim;
  Particles::ParticleHandler<dim, spacedim> particle_handler(tr,
                                                             mapping,
                                                             n_properties);

  std::vector<Point<dim>> particle_reference_locations =
    QGauss<dim>(3).get_points();

  Particles::Generators::regular_reference_locations(
    tr, particle_reference_locations, particle_handler, mapping);

  for (auto &particle : particle_handler)
    {
      particle.get_properties()[spacedim - 1] = particle.get_location()[0];
      deallog << "Before copying particle id " << particle.get_id()
              << " has first property " << particle.get_properties()[0]
              << " and last property "
              << particle.get_properties()[spacedim - 1] << " and position "
              << particle.get_location() << std::endl;
    }

  {
    Particles::ParticleHandler<dim, spacedim> particle_handler_copy;
    particle_handler_copy.copy_from(particle_handler);

    // Make sure the old particle handler and the property pool
    // are cleared. This catches problems if the new particles try to access
    // old memory addresses (this was a bug, fixed in
    // https://github.com/dealii/dealii/pull/11314)
    particle_handler.clear();

    for (const auto &particle : particle_handler_copy)
      deallog << "After copying particle id " << particle.get_id()
              << " has first property " << particle.get_properties()[0]
              << " and last property "
              << particle.get_properties()[spacedim - 1] << " and position "
              << particle.get_location() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

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
