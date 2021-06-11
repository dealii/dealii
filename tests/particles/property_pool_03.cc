// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// test the sorting of memory chunks inside the property pool

#include <deal.II/particles/property_pool.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"


void
test()
{
  {
    const int dim      = 2;
    const int spacedim = 2;

    const unsigned int                     n_properties = 3;
    Particles::PropertyPool<dim, spacedim> pool(n_properties);

    std::vector<
      std::vector<typename Particles::PropertyPool<dim, spacedim>::Handle>>
      particle_handles;
    particle_handles.resize(2);

    // Allocate some space in non-contigous locations
    particle_handles[1].push_back(pool.register_particle());
    particle_handles[0].push_back(pool.register_particle());
    particle_handles[1].push_back(pool.register_particle());
    particle_handles[0].push_back(pool.register_particle());
    particle_handles[1].push_back(pool.register_particle());
    particle_handles[0].push_back(pool.register_particle());

    unsigned int i = 0;
    for (auto &particles_in_cell : particle_handles)
      for (auto &particle : particles_in_cell)
        pool.set_id(particle, i++);

    for (const auto &particles_in_cell : particle_handles)
      for (const auto &particle : particles_in_cell)
        deallog << "Before removal unsorted ID: " << pool.get_id(particle)
                << ". Handle: " << particle << std::endl;

    pool.sort_memory_slots(particle_handles);

    for (const auto &particles_in_cell : particle_handles)
      for (const auto &particle : particles_in_cell)
        deallog << "Before removal sorted ID: " << pool.get_id(particle)
                << ". Handle: " << particle << std::endl;

    // Deallocate some space in the same way the particle handler would
    // remove particles
    pool.deregister_particle(particle_handles[1][1]);
    pool.deregister_particle(particle_handles[0][2]);

    particle_handles[1][1] = particle_handles[1].back();
    particle_handles[1].resize(particle_handles[1].size() - 1);
    particle_handles[0].resize(particle_handles[0].size() - 1);

    for (const auto &particles_in_cell : particle_handles)
      for (const auto &particle : particles_in_cell)
        deallog << "After removal unsorted ID: " << pool.get_id(particle)
                << ". Handle: " << particle << std::endl;

    pool.sort_memory_slots(particle_handles);

    for (const auto &particles_in_cell : particle_handles)
      for (const auto &particle : particles_in_cell)
        deallog << "After removal sorted ID: " << pool.get_id(particle)
                << ". Handle: " << particle << std::endl;

    for (auto &particles_in_cell : particle_handles)
      for (auto &particle : particles_in_cell)
        pool.deregister_particle(particle);
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test();
}
