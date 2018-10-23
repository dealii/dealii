// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#include <deal.II/base/incremental_function.h>

#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim, typename RangeNumberType>
  IncrementalFunction<dim, RangeNumberType>::IncrementalFunction(
    Function<dim, RangeNumberType> &base)
    : base(base)
    , delta_t(time_type())
    , values_old(base.n_components)
  {}



  template <int dim, typename RangeNumberType>
  void
  IncrementalFunction<dim, RangeNumberType>::set_decrement(
    const time_type delta)
  {
    delta_t = delta;
  }



  template <int dim, typename RangeNumberType>
  void
  IncrementalFunction<dim, RangeNumberType>::vector_value(
    const Point<dim> &       p,
    Vector<RangeNumberType> &values) const
  {
    // since we modify a mutable member variable, lock the
    // the data via a mutex
    std::lock_guard<std::mutex> lock(mutex);

    base.set_time(this->get_time());
    base.vector_value(p, values);

    base.set_time(this->get_time() - delta_t);
    base.vector_value(p, values_old);

    values -= values_old;
  }



  template <int dim, typename RangeNumberType>
  RangeNumberType
  IncrementalFunction<dim, RangeNumberType>::value(const Point<dim> &p,
                                                   unsigned int      comp) const
  {
    base.set_time(this->get_time());
    const RangeNumberType current = base.value(p, comp);

    base.set_time(this->get_time() - delta_t);
    const RangeNumberType old = base.value(p, comp);

    return current - old;
  }


// Explicit instantiations
#include "incremental_function.inst"
} // namespace Functions
DEAL_II_NAMESPACE_CLOSE
