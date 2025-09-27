// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Implement some serialization comparison specializations

#include <deal.II/differentiation/sd/symengine_optimizer.h>

#include "../../serialization/serialization.h"

using namespace dealii;
namespace SD = Differentiation::SD;

// Overload the comparison function in serialization.h
template <typename ReturnType>
bool
compare_optimizers(const SD::BatchOptimizer<ReturnType> &t1,
                   const SD::BatchOptimizer<ReturnType> &t2)
{
  return (t1.n_independent_variables() == t2.n_independent_variables()) &&
         (t1.n_dependent_variables() == t2.n_dependent_variables()) &&
         (t1.optimization_method() == t2.optimization_method()) &&
         (t1.optimization_flags() == t2.optimization_flags());
}

template <>
bool
compare(const SD::BatchOptimizer<float> &t1,
        const SD::BatchOptimizer<float> &t2)
{
  return compare_optimizers(t1, t2);
}

template <>
bool
compare(const SD::BatchOptimizer<double> &t1,
        const SD::BatchOptimizer<double> &t2)
{
  return compare_optimizers(t1, t2);
}


template <>
bool
compare(const SD::BatchOptimizer<std::complex<float>> &t1,
        const SD::BatchOptimizer<std::complex<float>> &t2)
{
  return compare_optimizers(t1, t2);
}

template <>
bool
compare(const SD::BatchOptimizer<std::complex<double>> &t1,
        const SD::BatchOptimizer<std::complex<double>> &t2)
{
  return compare_optimizers(t1, t2);
}



template <typename T>
void
verify_no_logging(const T &t1, T &t2)
{
  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << t1;
    // archive and stream closed when
    // destructors are called
  }

  // verify correctness of the
  // serialization
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);


    ia >> t2;

    AssertThrow(compare(t1, t2), ExcInternalError());
  }
}
