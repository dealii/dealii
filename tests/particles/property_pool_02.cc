// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
      deallog << ' ' << pool.get_properties(handle)[i];

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
