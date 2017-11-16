// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Like particle_02, but with particle properties.

#include "../tests.h"
#include <deal.II/particles/particle.h>
#include <deal.II/base/array_view.h>


template <int dim>
void test ()
{
  {
    const unsigned int n_properties_per_particle = 3;
    Particles::PropertyPool pool(n_properties_per_particle);

    Point<2> position;
    position(0) = 0.3;
    position(1) = 0.5;

    Point<2> reference_position;
    reference_position(0) = 0.2;
    reference_position(1) = 0.4;

    const types::particle_index index(7);

    std::vector<double> properties = {0.15,0.45,0.75};

    Particles::Particle<2> particle(position,reference_position,index);
    particle.set_property_pool(pool);
    particle.set_properties(ArrayView<double>(&properties[0],properties.size()));

    deallog << "Particle properties: "
            << std::vector<double>(particle.get_properties().begin(),particle.get_properties().end())
            << std::endl;

    const Particles::Particle<2> copy(particle);

    deallog << "Copy particle properties: "
            << std::vector<double>(copy.get_properties().begin(),copy.get_properties().end())
            << std::endl;

    deallog << "Old particle has properties before move: " << particle.has_properties() << std::endl;

    const Particles::Particle<2> moved_particle(std::move(particle));

    deallog << "Old particle has properties after move: " << particle.has_properties() << std::endl;

    deallog << "Moved particle properties: "
            << std::vector<double>(moved_particle.get_properties().begin(),moved_particle.get_properties().end())
            << std::endl;

  }

  deallog << "OK" << std::endl;
}



int main ()
{
  initlog();
  test<2>();
}
