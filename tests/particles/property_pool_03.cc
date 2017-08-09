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



// test a property pool that uses a custom type as storage container

#include "../tests.h"
#include <deal.II/particles/property_pool.templates.h>
#include <fstream>
#include <iomanip>



struct DataType
{
  double a[2];
  char c;
  double *pointer;
};

void test ()
{
  {
    double b[2] = {0.7,1.3};

    Particles::PropertyPool<DataType> pool;

    typename Particles::PropertyPool<DataType>::Handle handle = pool.allocate_properties_array();

    DataType &properties_in = pool.get_properties(handle)[0];
    properties_in.a[0] = 2.5;
    properties_in.a[1] = 2.0;
    properties_in.c = 'f';
    properties_in.pointer = &b[1];

    const DataType &properties_out = pool.get_properties(handle)[0];

    deallog << "Pool properties: " << properties_out.a[0] << ' '
            << properties_out.a[1] << ' '
            << properties_out.c << ' '
            << *properties_out.pointer
            << std::endl;

    pool.deallocate_properties_array(handle);
  }

  deallog << "OK" << std::endl;
}



int main ()
{
  initlog();
  ;

  test();
}
