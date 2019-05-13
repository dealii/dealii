// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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


// check that internal::VectorOperations::parallel_reduce works for start-end

#include <deal.II/base/index_set.h>

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

      for (unsigned int i = 0; i < size; ++i)
        val[i] = random_value<double>();


      internal::VectorOperations::MeanValue<Number> mean(val);

      Number sum_direct = 0.;
      internal::VectorOperations::parallel_reduce(
        mean, 0, size, sum_direct, thread_loop_partitioner);

      Number sum = 0.;
      // now break the size in chunks
      const unsigned int n_chunks   = 3;
      const unsigned int chunk_size = size / n_chunks;
      for (unsigned int i = 0; i <= n_chunks; ++i)
        {
          const unsigned int begin = i * chunk_size;
          const unsigned int end   = std::min((i + 1) * chunk_size, size);

          Number sum_i = 0.;
          internal::VectorOperations::parallel_reduce(
            mean, begin, end, sum_i, thread_loop_partitioner);
          sum += sum_i;
        }

      // check values:
      AssertThrow(std::fabs(sum - sum_direct) < 1e-6 * sum_direct,
                  ExcInternalError());

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
