// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
