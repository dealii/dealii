// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/incremental_function.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim, typename RangeNumberType>
  IncrementalFunction<dim, RangeNumberType>::IncrementalFunction(
    Function<dim, RangeNumberType> &base)
    : Function<dim, RangeNumberType>(base.n_components)
    , base(base)
    , delta_t(numbers::signaling_nan<time_type>())
    , values_old(base.n_components)
  {}



  template <int dim, typename RangeNumberType>
  void
  IncrementalFunction<dim, RangeNumberType>::set_decrement(
    const time_type delta)
  {
    Assert(delta >= 0.0,
           ExcMessage("The decrement must be set to a non-negative value."));
    delta_t = delta;
  }



  template <int dim, typename RangeNumberType>
  RangeNumberType
  IncrementalFunction<dim, RangeNumberType>::value(
    const Point<dim>  &p,
    const unsigned int comp) const
  {
    // since we modify a mutable member variable, lock the
    // the data via a mutex
    std::lock_guard<std::mutex> lock(mutex);

    // Cache the time state of the base class in case it has been changed
    // within the user code. We reset the wrapped function to the original
    // state once we're done with our own evaluations.
    const auto orig_time = base.get_time();

    base.set_time(this->get_time());
    const RangeNumberType current = base.value(p, comp);

    base.set_time(this->get_time() - delta_t);
    const RangeNumberType old = base.value(p, comp);

    // Reset wrapped function time setting
    base.set_time(orig_time);

    return current - old;
  }



  template <int dim, typename RangeNumberType>
  void
  IncrementalFunction<dim, RangeNumberType>::vector_value(
    const Point<dim>        &p,
    Vector<RangeNumberType> &values) const
  {
    // since we modify a mutable member variable, lock the
    // the data via a mutex
    std::lock_guard<std::mutex> lock(mutex);

    // Cache the time state of the base class in case it has been changed
    // within the user code. We reset the wrapped function to the original
    // state once we're done with our own evaluations.
    const auto orig_time = base.get_time();

    base.set_time(this->get_time());
    base.vector_value(p, values);

    base.set_time(this->get_time() - delta_t);
    base.vector_value(p, values_old);

    values -= values_old;

    // Reset wrapped function time setting
    base.set_time(orig_time);
  }


// Explicit instantiations
#include "base/incremental_function.inst"
} // namespace Functions
DEAL_II_NAMESPACE_CLOSE
