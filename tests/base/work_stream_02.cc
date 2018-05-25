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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
