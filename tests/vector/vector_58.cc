// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check Vector<double>::size()

#include <deal.II/lac/vector.h>

#include "../tests.h"

static unsigned int counter = 0;

template <typename Number>
void
fill(Vector<Number> &v)
{
  v = 0;
  for (unsigned int i = 0; i < v.size(); i++)
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
