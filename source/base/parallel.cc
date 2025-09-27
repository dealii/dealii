// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/exceptions.h>
#include <deal.II/base/parallel.h>

#ifdef DEAL_II_WITH_TASKFLOW
#  include <deal.II/base/multithread_info.h>

#  include <taskflow/algorithm/for_each.hpp>
#  include <taskflow/taskflow.hpp>
#endif

#ifdef DEAL_II_WITH_TBB
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/partitioner.h>
#else
#  include <boost/range/iterator_range.hpp>
#endif


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace VectorImplementation
  {
    // set minimum grain size. this value has been determined by experiments
    // with the actual dealii::Vector implementation and takes vectorization
    // that is done inside most functions of dealii::Vector into account
    // (without vectorization, we would land at around 1000). for smaller
    // values than this, the scheduling overhead will be significant. if the
    // value becomes too large, on the other hand, the work load cannot be
    // split into enough chunks and the problem becomes badly balanced. the
    // tests are documented at https://github.com/dealii/dealii/issues/2496
    unsigned int minimum_parallel_grain_size = 4096;
  } // namespace VectorImplementation


  namespace SparseMatrixImplementation
  {
    // set this value to 1/16 of the value of the minimum grain size of
    // vectors (a factor of 4 because we assume that sparse matrix-vector
    // products do not vectorize and the other factor of 4 because we expect
    // at least four entries per row). this rests on the fact that we have to
    // do a lot more work per row of a matrix than per element of a vector. it
    // could possibly be reduced even further but that doesn't appear worth it
    // any more for anything but very small matrices that we don't care that
    // much about anyway.
    unsigned int minimum_parallel_grain_size = 256;
  } // namespace SparseMatrixImplementation
} // namespace internal

namespace parallel
{
  namespace internal
  {
#ifdef DEAL_II_WITH_TBB
    TBBPartitioner::TBBPartitioner()
      : my_partitioner(std::make_shared<tbb::affinity_partitioner>())
      , in_use(false)
    {}



    TBBPartitioner::~TBBPartitioner()
    {
      AssertNothrow(in_use == false,
                    ExcInternalError(
                      "A vector partitioner goes out of scope, but "
                      "it appears to be still in use."));
    }



    std::shared_ptr<tbb::affinity_partitioner>
    TBBPartitioner::acquire_one_partitioner()
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (in_use)
        return std::make_shared<tbb::affinity_partitioner>();

      in_use = true;
      return my_partitioner;
    }



    void
    TBBPartitioner::release_one_partitioner(
      const std::shared_ptr<tbb::affinity_partitioner> &p)
    {
      if (p.get() == my_partitioner.get())
        {
          std::lock_guard<std::mutex> lock(mutex);
          in_use = false;
        }
    }
#else
    TBBPartitioner::TBBPartitioner() = default;
#endif
  } // namespace internal
} // namespace parallel



DEAL_II_NAMESPACE_CLOSE
