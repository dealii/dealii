// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test functions in namespace WorkStream

#include <deal.II/base/work_stream.h>

#include "../tests.h"


struct ScratchData
{};


struct CopyData
{
  unsigned int computed;
};


struct X
{
  void
  worker(const std::vector<unsigned int>::iterator &i,
         ScratchData &,
         CopyData &ad)
  {
    ad.computed = *i * 2;
  }

  void
  copier(const CopyData &ad)
  {
    deallog << ad.computed << std::endl;
  }
};


void
test()
{
  std::vector<unsigned int> v;
  for (unsigned int i = 0; i < 20; ++i)
    v.push_back(i);

  X x;
  WorkStream::run(
    v.begin(), v.end(), x, &X::worker, &X::copier, ScratchData(), CopyData());
}



int
main()
{
  initlog();

  test();
}
