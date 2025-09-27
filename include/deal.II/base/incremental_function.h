// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_incremental_function_h
#define dealii_incremental_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/mutex.h>

#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

namespace Functions
{
  /**
   * This class represents an incremental function. That is, given an arbitrary
   * function <code>f</code>, this class will return
   * <code>f(t) - f(t - delta_t)</code>, where <code>f(t)</code> denotes the
   * function evaluated at time <code>t</code> and, likewise, <code>f(t -
   * delta_t)</code> denotes the function evaluated at time <code>t -
   * delta_t</code>. The decrement <code>delta_t</code> is set by the method
   * set_decrement(). The main application of this class is to transform a given
   * Dirichlet boundary condition function into incremental form, as is
   * required by some implementations of non-linear solution schemes.
   *
   * @ingroup functions
   */
  template <int dim, typename RangeNumberType = double>
  class IncrementalFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * Export the value of the template parameter as a static member constant.
     * This is sometimes useful in the context of template programming.
     */
    static constexpr unsigned int dimension = dim;

    /**
     * The scalar-valued real type used for representing time.
     */
    using time_type = typename Function<dim, RangeNumberType>::time_type;

    /**
     * Constructor which wraps a given function @p base.
     *
     * @note This class stores a non-constant reference to @p base
     * and will call <code>base.set_time()</code> during evaluation
     * in order to evaluate the @p base class at any arbitrary time.
     * It is guaranteed that the temporal state of @p base is returned
     * to its original settings after each function evaluation in this
     * class.
     */
    IncrementalFunction(Function<dim, RangeNumberType> &base);

    /**
     * Return the value of the function at the given point.
     *
     * Unless there is only one component (i.e. the function is scalar), you
     * should state the component you want to have evaluated. By default, the
     * value of the first component is computed.
     */
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return all components of a vector-valued function at a given point.
     *
     * It is required that the @p values vector have the correct size before
     * this function is called.
     */
    virtual void
    vector_value(const Point<dim>        &p,
                 Vector<RangeNumberType> &values) const override;

    /**
     * Set the time decrement.
     *
     * It is expected that this value be positive.
     */
    void
    set_decrement(const time_type delta_t);

  private:
    /**
     * A reference to the function being wrapped.
     */
    Function<dim, RangeNumberType> &base;

    /**
     * The time decrement.
     */
    time_type delta_t;

    /**
     * An auxiliary vector to store values.
     */
    mutable Vector<RangeNumberType> values_old;

    /**
     * Thread mutex for supporting evaluation in multi-threaded contexts.
     */
    mutable Threads::Mutex mutex;
  };

} // namespace Functions


DEAL_II_NAMESPACE_CLOSE

#endif
