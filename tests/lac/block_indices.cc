// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


// Test if block indices are handled properly

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_indices.h>

#include <fstream>


void test (const BlockIndices &idx)
{
  const unsigned int n = idx.size();
  deallog << "sizes: " << idx << std::endl;
  deallog << "start:   ";
  for (unsigned i=0; i<n; ++i)
    deallog << ' ' << idx.block_start(i);

  deallog << std::endl;

  for (unsigned int i=0; i<idx.total_size(); ++i)
    {
      const unsigned int b = idx.global_to_local(i).first;
      const unsigned int j = idx.global_to_local(i).second;
      deallog << ' '<< i << ':' << b << ':' << j;
    }

  deallog << std::endl;

  for (unsigned int b=0; b<n; ++b)
    for (unsigned int j=0; j<idx.block_size(b); ++j)
      {
        const unsigned int i = idx.local_to_global(b,j);
        deallog << ' '<< i << ':' << b << ':' << j;
      }

  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  BlockIndices bi0;
  deallog << "empty: " << bi0 << std::endl;
  bi0.push_back(3);
  deallog << "push:  " << bi0 << std::endl;
  bi0.push_back(2);
  deallog << "push:  " << bi0 << std::endl;
  bi0.reinit(0,0);
  deallog << "empty: " << bi0 << std::endl;

  BlockIndices bi1(3);
  test(bi1);
  bi1.reinit(3,4);
  test(bi1);
  bi1.push_back(2);
  deallog << "push: " << bi1 << std::endl;
  bi1.push_back(5);
  deallog << "push: " << bi1 << std::endl;
  bi1.push_back(4);
  test(bi1);

  std::vector<types::global_dof_index> v(4);
  for (unsigned int i=0; i<v.size(); ++i)
    v[i] = 4-i;

  BlockIndices bi2(v);
  test(bi2);
}
