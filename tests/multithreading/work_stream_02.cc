// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test functions in namespace WorkStream

#include <deal.II/base/work_stream.h>

#include "../tests.h"


struct ScratchData
{};


struct CopyData
{
  unsigned int computed;
};


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


void
test()
{
  std::vector<unsigned int> v;
  for (unsigned int i = 0; i < 20; ++i)
    v.push_back(i);

  WorkStream::run(
    v.begin(), v.end(), &worker, &copier, ScratchData(), CopyData());
}



int
main()
{
  initlog();

  test();
}
