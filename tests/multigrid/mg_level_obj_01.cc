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


// check MGLevelObject::apply()

#include "../tests.h"
#include <deal.II/base/mg_level_object.h>

#include <algorithm>

using namespace std;


template <class T>
void
check(MGLevelObject<T> &obj)
{
  obj.apply([&] (const unsigned int lvl, T & value)
  {
    value = (T)lvl;
  });

  obj.apply([&] (const unsigned int lvl, T & value)
  {
    deallog << "lvl: " << lvl << " value: " << value << std::endl;
  });
}

int
main()
{
  initlog();

  MGLevelObject<double> o(2,4);
  check(o);
  o.resize(0,1);
  check(o);
}
