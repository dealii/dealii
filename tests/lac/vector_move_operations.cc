// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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


#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();

  {
    Vector<double> temp(10);
    Vector<double> w(std::move(temp));
    deallog << "move constructor:  " << w << std::endl;
  }

  deallog << std::endl;

  {
    Vector<double> u(10);
    for (unsigned int i = 0; i < u.size(); ++i)
      u[i] = (double)(i + 1);

    deallog << "vector:          " << u << std::endl;

    Vector<double> v;
    v = u;
    deallog << "copy assignment: " << v << std::endl;
    deallog << "old object:      " << u << std::endl;

    v = std::move(u);
    deallog << "move assignment: " << v << std::endl;
    deallog << "old object size: " << u.size() << std::endl;

    // and swap again with different sizes
    u.reinit(0);
    u = std::move(v);
    deallog << "move assignment: " << u << std::endl;
    deallog << "old object size: " << v.size() << std::endl;
  }
}
