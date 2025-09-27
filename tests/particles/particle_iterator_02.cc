// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like particle_iterator_01, but tests state of the iterator.

#include <deal.II/base/array_view.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include "../tests.h"


template <int dim>
void
test()
{
  {
    Triangulation<dim> tr;

    GridGenerator::hyper_cube(tr);
    MappingQ<dim> mapping(1);

    const unsigned int           n_properties_per_particle = 3;
    Particles::PropertyPool<dim> pool(n_properties_per_particle);

    Point<dim> position;
    position[0] = 0.3;

    if (dim > 1)
      position[1] = 0.5;

    Point<dim> reference_position;
    reference_position[0] = 0.2;

    if (dim > 1)
      reference_position[1] = 0.4;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15, 0.45, 0.75};

    typename Particles::ParticleHandler<dim>::particle_container
      particle_container;

    auto particle = pool.register_particle();
    particle_container.emplace_back(
      std::vector<typename Particles::PropertyPool<dim>::Handle>(1, particle),
      tr.begin_active());
    particle_container.emplace_front(
      std::vector<typename Particles::PropertyPool<dim>::Handle>(), tr.end());
    particle_container.emplace_back(
      std::vector<typename Particles::PropertyPool<dim>::Handle>(), tr.end());

    pool.set_location(particle, position);
    pool.set_reference_location(particle, reference_position);
    pool.set_id(particle, index);
    auto particle_properties = pool.get_properties(particle);

    std::copy(properties.begin(),
              properties.end(),
              particle_properties.begin());

    Particles::ParticleIterator<dim> particle_begin(
      ++particle_container.begin(), pool, 0);
    Particles::ParticleIterator<dim> particle_end(--particle_container.end(),
                                                  pool,
                                                  0);
    Particles::ParticleIterator<dim> particle_nonexistent1(
      ++particle_container.begin(), pool, 1);
    Particles::ParticleIterator<dim> particle_nonexistent2(
      --particle_container.end(), pool, 1);
    Particles::ParticleIterator<dim> particle_invalid;

    Assert(particle_begin->state() == IteratorState::valid, ExcInternalError());
    Assert(particle_end->state() == IteratorState::past_the_end,
           ExcInternalError());
    Assert(particle_nonexistent1->state() == IteratorState::invalid,
           ExcInternalError());
    Assert(particle_nonexistent2->state() == IteratorState::invalid,
           ExcInternalError());
    Assert(particle_invalid->state() == IteratorState::invalid,
           ExcInternalError());

    Particles::ParticleIterator<dim> particle_iterator = particle_begin;
    ++particle_iterator;

    Assert(particle_iterator->state() == IteratorState::past_the_end,
           ExcInternalError());

    particle_iterator = particle_begin;
    --particle_iterator;

    Assert(particle_iterator->state() == IteratorState::past_the_end,
           ExcInternalError());

    for (particle_iterator = particle_begin; particle_iterator != particle_end;
         ++particle_iterator)
      {
        deallog << "Particle position: " << particle_iterator->get_location()
                << ". Iterator valid: "
                << (particle_iterator->state() == IteratorState::valid)
                << std::endl;
      }

    Assert(particle_iterator->state() == IteratorState::past_the_end,
           ExcInternalError());

    pool.deregister_particle(particle);
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
  test<3>();
}
