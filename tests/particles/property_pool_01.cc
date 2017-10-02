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



// check the creation, simplest usage, and destruction of a property pool

#include "../tests.h"
#include <deal.II/particles/property_pool.h>
#include <fstream>
#include <iomanip>


void test ()
{
  {
    Particles::PropertyPool pool;

    typename Particles::PropertyPool::Handle handle = pool.allocate_properties_array();

    pool.get_properties(handle)[0] = 2.5;

    deallog << "Pool properties: " << pool.get_properties(handle)[0] << std::endl;

    pool.deallocate_properties_array(handle);
  }

  deallog << "OK" << std::endl;
}



int main ()
{
  initlog();
  test();
}
