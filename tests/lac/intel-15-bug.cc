// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

// intel 15.0 compiler bug. This is a simplification of tests/lac/vector-vector
// It only triggers when using TBB (it will use two threads) and only with
// SIMD for long double. We now use dealii::parallel::internal::EnableOpenMPSimdFor
// so the test passes.

#include <deal.II/base/parallel.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

typedef int size_type;

template <typename Number>
struct Vectorization_add_v
{
  Number *val;
  Number *v_val;

  void
  operator()(const tbb::blocked_range<size_type> &range) const
  {
    if (dealii::parallel::internal::EnableOpenMPSimdFor<Number>::value)
      {
        DEAL_II_OPENMP_SIMD_PRAGMA
        for (size_type i = range.begin(); i < range.end(); ++i)
          val[i] += v_val[i];
      }
    else
      {
        for (size_type i = range.begin(); i < range.end(); ++i)
          val[i] += v_val[i];
      }
  }
};

const unsigned int N = 3;

void
check()
{
  std::vector<long double> d1(N), d2(N);
  for (unsigned int i = 0; i < N; ++i)
    {
      d1[i] = 1.0;
      d2[i] = i;
    };

  Vectorization_add_v<long double> vector_add;
  vector_add.val   = &d1[0];
  vector_add.v_val = &d2[0];

  if (1)
    {
      //fails:
      tbb::parallel_for(tbb::blocked_range<size_type>(0, N, 2),
                        vector_add,
                        tbb::auto_partitioner());
    }
  else
    {
      // works:
      vector_add(tbb::blocked_range<int>(0, 1, 2));
      vector_add(tbb::blocked_range<int>(1, 3, 2));
    }

  deallog << "result: " << d1[N - 1] << " should be 3" << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check();
}
