// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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


// test WorkStream with empty functions

#include <deal.II/base/work_stream.h>

#include "../tests.h"


struct ScratchData
{};


void
foo(const std::vector<unsigned int>::iterator, ScratchData &, unsigned int &)
{}

void
bar(const unsigned int &)
{}

void
test()
{
  std::vector<unsigned int> v;
  for (unsigned int i = 0; i < 20; ++i)
    v.push_back(i);

  // first run with only a worker
  WorkStream::run(v.begin(),
                  v.end(),
                  &foo,
                  std::function<void(const unsigned int &)>(),
                  ScratchData(),
                  0U);

  // next run with only a copier
  WorkStream::run(v.begin(),
                  v.end(),
                  std::function<void(const std::vector<unsigned int>::iterator,
                                     ScratchData &,
                                     unsigned int &)>(),
                  &bar,
                  ScratchData(),
                  0U);
}



int
main()
{
  initlog();

  test();
}
