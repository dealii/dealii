// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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


// test functions in namespace WorkStream
// This test is based off base/work_stream_01.cc, but for iterator ranges

#include <deal.II/base/work_stream.h>

#include <deal.II/grid/filtered_iterator.h>

#include "../tests.h"


struct ScratchData
{};


struct CopyData
{
  unsigned int computed;
};


using IteratorType      = typename std::vector<unsigned int>::iterator;
using IteratorRangeType = IteratorRange<IteratorType>;

static_assert(
  std::is_same<typename IteratorRangeType::iterator, IteratorType>::value,
  "Iterator types not the same");

struct X
{
  void
  worker(const IteratorType &i, ScratchData &, CopyData &ad)
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

  const IteratorRangeType iterator_range(v.begin(), v.end());

  X x;
  WorkStream::run(
    iterator_range, x, &X::worker, &X::copier, ScratchData(), CopyData());
}



int
main()
{
  initlog();

  test();
}
