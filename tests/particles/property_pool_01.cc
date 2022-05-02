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



// check the creation, simplest usage, and destruction of a property pool

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

    Particles::PropertyPool<dim, spacedim> pool(1);

    typename Particles::PropertyPool<dim, spacedim>::Handle handle =
      pool.register_particle();

    pool.get_properties(handle)[0] = 2.5;

    deallog << "Pool properties: " << pool.get_properties(handle)[0]
            << std::endl;

    pool.deregister_particle(handle);
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test();
}
