// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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


#include <deal.II/base/parallel.h>


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
      std::shared_ptr<tbb::affinity_partitioner> &p)
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
