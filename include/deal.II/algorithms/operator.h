// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


#ifndef dealii_operator_h
#define dealii_operator_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>

#include <deal.II/base/event.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace containing numerical algorithms in a unified form.
 *
 * All algorithmic classes in this namespace are derived from either Operator
 * or OutputOperator, depending on whether they return a value or not. See the
 * documentation of those classes for more detailed information on how to use
 * them.
 *
 * @author Guido Kanschat
 * @date 2012, 2013
 */
namespace Algorithms
{
  /**
   * @todo Update this documentation and the one of Operator
   *
   * The abstract base class of all algorithms in this library. An operator is
   * an object with an operator(), which transforms a set of named vectors
   * into another set of named vectors.
   *
   * Furthermore, an operator can be notified of parameter changes by the
   * calling routine. The outer iteration can notify() the Operator of an
   * Event, which could be for instance a change of mesh, a different time
   * step size or too slow convergence of Newton's method, which would then
   * trigger reassembling of a matrix or similar things.
   *
   * <h3>Usage for nested iterations</h3>
   *
   * This is probably the most prominent use for Operator, where an outer
   * iterative method calls an inner solver and so on. Typically, the
   * innermost method in such a nested system will have to compute a residual
   * using values from all outer iterations. Since the depth and order of such
   * a nesting is hardly predictable when designing a general tool, we use
   * AnyData to access these vectors. Typically, the first vector in
   * <tt>out</tt> contains the start vector when operator()() is called, and
   * the solution when the function returns. The object <tt>in</tt> is
   * providing additional information and forwarded to the inner Operator
   * objects of the nested iteration.
   *
   * @author Guido Kanschat
   * @date 2014
   */
  class OperatorBase : public Subscriptor
  {
  public:
    /**
     * The virtual destructor.
     */
    virtual ~OperatorBase() override = default;

    /**
     * The actual operation, which is implemented in a derived class.
     */
    virtual void
    operator()(AnyData &out, const AnyData &in) = 0;

    /**
     * Register an event triggered by an outer iteration.
     */
    virtual void
    notify(const Event &);

    /**
     * Clear all #notifications.
     */
    void
    clear_events();

  protected:
    /**
     * Accumulate events here. If any of those is set, the function solve() of
     * a terminal application must take care of reassembling the matrix.
     */
    Event notifications;
  };

  /**
   * An unary operator base class, intended to output the vectors in AnyData
   * in each step of an iteration.
   *
   * @author Guido Kanschat, 2010
   */
  template <typename VectorType>
  class OutputOperator : public Subscriptor
  {
  public:
    /**
     * Constructor initializing member variables with invalid data.
     */
    OutputOperator();

    /**
     * The copy constructor is deleted since objects of this class
     * should not be copied.
     */
    OutputOperator(const OutputOperator<VectorType> &) = delete;

    /**
     * Empty virtual destructor.
     */
    virtual ~OutputOperator() override = default;

    /**
     * Set the stream @p os to which data is written. If no stream is selected
     * with this function, data goes to @p deallog.
     */
    void
    initialize_stream(std::ostream &stream);

    /**
     * Set the current step.
     */
    void
    set_step(const unsigned int step);

    /**
     * Output all the vectors in AnyData.
     */
    virtual OutputOperator<VectorType> &
    operator<<(const AnyData &vectors);

  protected:
    unsigned int step;

  private:
    std::ostream *os;
  };

  template <typename VectorType>
  inline void
  OutputOperator<VectorType>::set_step(const unsigned int s)
  {
    step = s;
  }


  /**
   * Set the step number in OutputOperator by shifting an integer value.
   *
   * @relatesalso OutputOperator
   */
  template <typename VectorType>
  inline OutputOperator<VectorType> &
  operator<<(OutputOperator<VectorType> &out, unsigned int step)
  {
    out.set_step(step);
    return out;
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
