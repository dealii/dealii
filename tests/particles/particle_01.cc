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



// check the creation and destruction of particles

#include "../tests.h"
#include <deal.II/particles/particle.h>


template <int dim>
void test ()
{
  {
    Particles::Particle<dim> particle;

    deallog << "Particle location: " << particle.get_location() << std::endl;

    Point<dim> position;
    position(0) = 1.0;
    particle.set_location(position);

    deallog << "Particle location: " << particle.get_location() << std::endl;
  }

  deallog << "OK" << std::endl;
}



int main ()
{
  initlog();
  test<2>();
  test<3>();
}
