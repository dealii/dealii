// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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



// check the particle generation using the
// Particles::Generator::regular_reference_locations function using a nonlinear
// mapping in a hyper shell. The particles are generated at Gauss quadrature
// points.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  {
    Triangulation<dim, spacedim> tr;

    GridGenerator::hyper_shell(tr, Point<dim>(), 0.5, 1.0);

    MappingQGeneric<dim, spacedim> mapping(4);

    Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

    std::vector<Point<dim>> particle_reference_locations =
      QGauss<dim>(3).get_points();

    Particles::Generators::regular_reference_locations(
      tr, particle_reference_locations, particle_handler, mapping);

    deallog << "Particle number: " << particle_handler.n_global_particles()
            << std::endl;

    for (const auto &particle : particle_handler)
      {
        deallog << "Particle location: " << particle.get_location()
                << std::endl;
        deallog << "Particle reference location: "
                << particle.get_reference_location() << std::endl;
      }
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
