// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
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
  std::is_same_v<typename IteratorRangeType::iterator, IteratorType>,
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
