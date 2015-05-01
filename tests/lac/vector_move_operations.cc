// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>

#define PRINTME(name, var) \
  deallog << "Block vector: " name << ":" << std::endl; \
  for (unsigned int i = 0; i < var.n_blocks(); ++i) \
    deallog << "[block " << i << " ]  " << var.block(i);

int main()
{
  initlog();

  {
    Vector<double> w(Vector<double>(10));
    deallog << "move constructor:  " << w;
  }

  deallog << std::endl;

  {
    Vector<double> u(10);
    for (unsigned int i = 0; i < u.size(); ++i)
      u[i] = (double)(i+1);

    deallog << "vector:          " << u;

    Vector<double> v;
    v = u;
    deallog << "copy assignment: " << v;
    deallog << "old object:      " << u;

    v = std::move(u);
    deallog << "move assignment: " << v;
    deallog << "old object size: " << u.size() << std::endl;

    // and swap again with different sizes
    u = std::move(v);
    deallog << "move assignment: " << u;
    deallog << "old object size: " << v.size() << std::endl;
  }

  deallog << std::endl;

  {
    BlockVector<double> w(BlockVector<double>(5, 2));
    PRINTME("move constructor", w);
  }

  deallog << std::endl;

  {
    BlockVector<double> u(5, 2);
    for (unsigned int i = 0; i < 5; ++i)
      for (unsigned int j = 0; j < 2; ++j)
        u.block(i)[j] = (double)(10 * i + j);

    PRINTME("BlockVector", u);

    BlockVector<double> v;
    v = u;
    PRINTME("copy assignemnt", v);
    PRINTME("old object", u);

    v = std::move(u);
    PRINTME("move assignemnt", v);
    deallog << "old object size: " << u.n_blocks() << std::endl;

    // and swap again with different sizes
    u = std::move(v);
    PRINTME("move assignemnt", u);
    deallog << "old object size: " << v.n_blocks() << std::endl;
  }
}


