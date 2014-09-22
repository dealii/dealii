// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef __deal2__operator_h
#define __deal2__operator_h

#include <deal.II/base/config.h>
#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/named_data.h>
#include <deal.II/base/event.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace containing numerical algorithms in a unified form.
 *
 * All algorithmic classes in this namespace are derived from either
 * Operator or OutputOperator, depending on whether they return a
 * value or not. See the documentation of those classes for more
 * detailed information on how to use them.
 *
 * @author Guido Kanschat
 * @date 2012, 2013
 */
namespace Algorithms
{
  /**
   * @todo Update this documentation and the one of Operator
   *
   * The abstract base class of all algorithms in this library. An
   * operator is an object with an operator(), which transforms a set
   * of named vectors into another set of named vectors.
   *
   * Furthermore, an operator can be notified of parameter changes by
   * the calling routine. The outer iteration can notify() the Operator
   * of an Event, which could be for instance a change of mesh, a
   * different time step size or too slow convergence of Newton's
   * method, which would then trigger reassembling of a matrix or
   * similar things.
   *
   * <h3>Usage for nested iterations</h3>
   *
   * This is probably the most prominent use for Operator, where an
   * outer iterative method calls an inner solver and so on. Typically,
   * the innermost method in such a nested system will have to compute a
   * residual using values from all outer iterations. Since the depth
   * and order of such a nesting is hardly predictable when designing a
   * general tool, we use AnyData to access these vectors. Typically,
   * the first vector in <tt>out</tt> contains the start vector when
   * operator()() is called, and the solution when the function
   * returns. The object <tt>in</tt> is providing additional information
   * and forwarded to the inner Operator objects of the nested
   * iteration.
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
    ~OperatorBase();

    /**
     * The actual operation, which
     * is implemented in a derived class.
     */
    virtual void operator() (AnyData &out, const AnyData &in) = 0;

    /**
     * Register an event triggered
     * by an outer iteration.
     */
    virtual void notify(const Event &);
    /**
     * Clear all #notifications.
     */
    void clear_events();
  protected:
    /**
     * Accumulate events here. If any of
     * those is set, the function
     * solve() of a terminal
     * application must take care
     * of reassembling the matrix.
     */
    Event notifications;

  };

  /**
   * @deprecated This class has been replaced by OperatorBase.
   *
   * The abstract base class of all algorithms in this library. An
   * operator is an object with an operator(), which transforms a set
   * of named vectors into another set of named vectors.
   *
   * Furthermore, an operator can be notified of parameter changes by
   * the calling routine. The outer iteration can notify() the Operator
   * of an Event, which could be for instance a change of mesh, a
   * different time step size or too slow convergence of Newton's
   * method, which would then trigger reassembling of a matrix or
   * similar things.
   *
   * <h3>Usage for nested iterations</h3>
   *
   * This is probably the most prominent use for Operator, where an
   * outer iterative method calls an inner solver and so on. Typically,
   * the innermost method in such a nested system will have to compute a
   * residual using values from all outer iterations. Since the depth
   * and order of such a nesting is hardly predictable when designing a
   * general tool, we use NamedData to access these vectors. Typically,
   * the first vector in <tt>out</tt> contains the start vector when
   * operator()() is called, and the solution when the function
   * returns. The object <tt>in</tt> is providing additional information
   * and forwarded to the inner Operator objects of the nested
   * iteration.
   *
   * @author Guido Kanschat, 2010
   */
  template <class VECTOR>
  class Operator : public OperatorBase
  {
  public:
    Operator();

    /**
     * Implementation of the function in the base class in order to do
     * compatibility conversions between the old and the new
     * interface.
     */
    virtual void operator() (AnyData &out, const AnyData &in);

    /**
     * @deprecated It is in particular this function which should not be used anymore.
     *
     * The actual operation, which
     * is implemented in a derived class.
     */
    virtual void operator() (NamedData<VECTOR *> &out, const NamedData<VECTOR *> &in);

    /**
     * Set this true to avoid compatibility warnings.
     */
    bool silent_compatibility;

  private:
    /**
     * While we are providing compatibility functions to the old
     * interface, this variable will ensure there is no endless loop.
     */
    bool compatibility_flag;
  };

  /**
   * An unary operator base class, intended to output the vectors in
   * NamedData in each step of an iteration.
   *
   * @author Guido Kanschat, 2010
   */
  template <class VECTOR>
  class OutputOperator : public Subscriptor
  {
    OutputOperator(const OutputOperator<VECTOR> &);
  public:
    OutputOperator ();
    /**
     * Empty virtual destructor.
     */
    virtual ~OutputOperator();

    /**
     * Set the stream @p os to
     * which data is written. If
     * no stream is selected with
     * this function, data goes
     * to @p deallog.
     */
    void initialize_stream(std::ostream &stream);
    /**
     * Set the current step.
     */
    OutputOperator<VECTOR> &operator<< (unsigned int step);

    /**
     * Output all the vectors in AnyData.
     */
    virtual OutputOperator<VECTOR> &operator<< (const AnyData &vectors);

    /**
     * @deprecated Output all the vectors in NamedData.
     */
    OutputOperator<VECTOR> &operator<< (const NamedData<VECTOR *> &vectors);
  protected:
    unsigned int step;
  private:
    std::ostream *os;
  };

  template <class VECTOR>
  OutputOperator<VECTOR> &
  OutputOperator<VECTOR>::operator<< (unsigned int s)
  {
    step = s;
    return *this;
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
