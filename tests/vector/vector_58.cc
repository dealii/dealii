// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Vector<double>::size()

#include <deal.II/lac/vector.h>

#include "../tests.h"

static unsigned int counter = 0;

template <typename Number>
void
fill(Vector<Number> &v)
{
  v = 0;
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = counter + i * 2;

  ++counter;
}


template <typename Number>
void
test(const std::vector<unsigned int> &size_sequence)
{
  Vector<Number> v, v_old;
  for (unsigned int j = 0; j < size_sequence.size(); ++j)
    {
      const unsigned int s = size_sequence[j];
      if (v.size() == 0)
        {
          v.reinit(s);
        }
      else
        {
          const unsigned int check_s = (s > v.size() ? v.size() : s);
          v.grow_or_shrink(s);
          for (unsigned int i = 0; i < check_s; ++i)
            AssertThrow(v(i) == v_old(i),
                        ExcMessage("s=" + std::to_string(s) +
                                   " i=" + std::to_string(i) + " " +
                                   std::to_string(v(i)) +
                                   "!=" + std::to_string(v_old(i))));

          for (unsigned int i = check_s; i < s; ++i)
            AssertThrow(v(i) == 0.,
                        ExcMessage("s=" + std::to_string(s) +
                                   " i=" + std::to_string(i) + " " +
                                   std::to_string(v(i)) + "!=0"));
        }

      fill(v);
      v_old.reinit(s);
      v_old = v;
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  const std::vector<unsigned int> size_sequence = {
    {1,   2,  3,  5,  4,  7,  9,  7,   11,  15,
     220, 19, 18, 17, 16, 40, 35, 129, 300, 287}};
  test<double>(size_sequence);
}
