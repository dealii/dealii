// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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

    std::vector<typename Particles::PropertyPool<dim, spacedim>::Handle>
      particle_handles;

    // Allocate some space in non-contigous locations
    particle_handles.push_back(pool.register_particle());
    particle_handles.insert(particle_handles.begin(), pool.register_particle());
    particle_handles.push_back(pool.register_particle());
    particle_handles.insert(particle_handles.begin(), pool.register_particle());
    particle_handles.push_back(pool.register_particle());
    particle_handles.insert(particle_handles.begin(), pool.register_particle());

    unsigned int i = 0;
    for (auto &particle : particle_handles)
      pool.set_id(particle, i++);

    for (const auto &particle : particle_handles)
      deallog << "Before removal unsorted ID: " << pool.get_id(particle)
              << ". Handle: " << particle << std::endl;

    pool.sort_memory_slots(particle_handles);

    typename Particles::PropertyPool<dim, spacedim>::Handle h = 0;
    for (auto &particle : particle_handles)
      particle = h++;

    for (const auto &particle : particle_handles)
      deallog << "Before removal sorted ID: " << pool.get_id(particle)
              << ". Handle: " << particle << std::endl;

    // Deallocate some space in the same way the particle handler would
    // remove particles
    pool.deregister_particle(particle_handles[2]);
    pool.deregister_particle(particle_handles[3]);

    particle_handles[2] = particle_handles.back();
    particle_handles.resize(particle_handles.size() - 1);

    particle_handles[3] = particle_handles.back();
    particle_handles.resize(particle_handles.size() - 1);

    for (const auto &particle : particle_handles)
      deallog << "After removal unsorted ID: " << pool.get_id(particle)
              << ". Handle: " << particle << std::endl;

    pool.sort_memory_slots(particle_handles);

    h = 0;
    for (auto &particle : particle_handles)
      particle = h++;

    for (const auto &particle : particle_handles)
      deallog << "After removal sorted ID: " << pool.get_id(particle)
              << ". Handle: " << particle << std::endl;

    for (auto &particle : particle_handles)
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
