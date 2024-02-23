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


// check that internal::VectorOperations::parallel_for works for start-end

#include <deal.II/base/index_set.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/vector_operations_internal.h>

#include "../tests.h"



template <typename Number>
void
check()
{
  for (unsigned int test = 0; test < 5; ++test)
    {
      const unsigned int size = 17 + test * 1101;

      std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        thread_loop_partitioner;
      thread_loop_partitioner.reset(
        new ::dealii::parallel::internal::TBBPartitioner());

      Number *val;
      Utilities::System::posix_memalign((void **)&val,
                                        64,
                                        sizeof(Number) * size);

      const Number                                   s = numbers::PI;
      internal::VectorOperations::Vector_set<Number> setter(s, val);

      // now break the size in chunks
      const unsigned int n_chunks   = 3;
      const unsigned int chunk_size = size / n_chunks;
      for (unsigned int i = 0; i <= n_chunks; ++i)
        {
          const unsigned int begin = i * chunk_size;
          const unsigned int end   = std::min((i + 1) * chunk_size, size);
          internal::VectorOperations::parallel_for(setter,
                                                   begin,
                                                   end,
                                                   thread_loop_partitioner);
        }

      // check values:
      for (unsigned int i = 0; i < size; ++i)
        AssertThrow(val[i] == s, ExcInternalError());

      free(val);
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check<float>();
  check<double>();
  check<long double>();
  deallog << "OK" << std::endl;
}
