// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// test a property pool that allocates more than one property per chunk

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

    typename Particles::PropertyPool<dim, spacedim>::Handle handle =
      pool.register_particle();

    pool.get_properties(handle)[0] = 1.2;
    pool.get_properties(handle)[1] = 2.5;
    pool.get_properties(handle)[2] = 2.7;


    deallog << "Pool properties:";

    for (unsigned int i = 0; i < pool.get_properties(handle).size(); ++i)
      deallog << " " << pool.get_properties(handle)[i];

    deallog << std::endl;

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
