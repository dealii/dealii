// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



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
