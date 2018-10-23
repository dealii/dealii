// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_incremental_function_h
#define dealii_incremental_function_h


#include <deal.II/base/function.h>
#include <deal.II/base/thread_management.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class Vector;

namespace Functions
{
  /**
   * This class represents an incremental function. That is, given arbitrary
   * function <code>func</code>, this class will return
   * <code>f(t)-f(t-delta)</code>. The decrement is set by the method
   * set_decrement(). The main application is to transform a given Dirichlet
   * boundary condition function into the incremental form.
   *
   * @ingroup functions
   * @author Denis Davydov, Jean-Paul Pelteret, 2018
   */
  template <int dim, typename RangeNumberType = double>
  class IncrementalFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * Export the value of the template parameter as a static member constant.
     * Sometimes useful for some expression template programming.
     */
    static const unsigned int dimension = dim;

    /**
     * The scalar-valued real type used for representing time.
     */
    using time_type = typename Function<dim, RangeNumberType>::time_type;

    /**
     * Constructor which wraps a given function @p base.
     *
     * @note this class stores a non-constant reference to @p base
     * and will call <code>base.set_time()</code> during evaluation.
     */
    IncrementalFunction(Function<dim, RangeNumberType> &base);

    /**
     * Virtual destructor
     */
    virtual ~IncrementalFunction() = default;

    /**
     * Return the value of the function at the given point. Unless there is only
     * one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the first
     * component.
     */
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Return all components of a vector-valued function at a given point.
     *
     * <tt>values</tt> shall have the right size beforehand.
     */
    virtual void
    vector_value(const Point<dim> &       p,
                 Vector<RangeNumberType> &values) const override;

    /**
     * Set (positive) time decrement.
     */
    void
    set_decrement(const time_type delta_t);

  private:
    /**
     * Reference to the function being wrapped.
     */
    Function<dim, RangeNumberType> &base;

    /**
     * Time decrement.
     */
    time_type delta_t;

    /**
     * Auxiliary vector to store values
     */
    mutable Vector<RangeNumberType> values_old;

    /**
     * Thread mutex
     */
    mutable Threads::Mutex mutex;
  };

} // namespace Functions


DEAL_II_NAMESPACE_CLOSE

#endif
