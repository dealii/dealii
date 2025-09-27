// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

    // Allocate some space in non-contiguous locations
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
