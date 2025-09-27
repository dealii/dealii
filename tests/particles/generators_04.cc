// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the creation and destruction of a random particle in a cell

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

    GridGenerator::hyper_cube(tr);
    MappingQ<dim, spacedim> mapping(1);

    std::mt19937 random_number_generator(time(nullptr));

    const Particles::Particle<dim, spacedim> particle =
      Particles::Generators::random_particle_in_cell(tr.begin_active(),
                                                     42,
                                                     random_number_generator,
                                                     mapping);

    const Point<dim> p_unit =
      mapping.transform_real_to_unit_cell(tr.begin_active(),
                                          particle.get_location());

    deallog << "Particle is inside cell 1: "
            << GeometryInfo<dim>::is_inside_unit_cell(p_unit) << std::endl;
    deallog << "Particle id: " << particle.get_id() << std::endl;
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
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
