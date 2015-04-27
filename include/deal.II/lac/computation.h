// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

#ifndef __deal2__computation_h
#define __deal2__computation_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector_memory.h>

#ifdef DEAL_II_WITH_CXX11

#include <functional>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:

template <typename Range, typename Domain> class LinearOperator;
template <typename Range = Vector<double> > class Computation;


/**
 * A class to store a computation.
 *
 * The class essentially consists of <code>std::function</code> objects
 * that store the knowledge of how to generate the result of a computation
 * and store it in a vector:
 * @code
 *   std::function<void(Range &)> result;
 *   std::function<void(Range &)> result_add;
 * @endcode
 *
 * Similar to the LinearOperator class it also has knowledge about how to
 * initialize a vector of the @p Range space:
 * @code
 *   std::function<void(Range &, bool)> reinit_vector;
 * @endcode
 *
 * The primary purpose of this class is to allow lazy evaluation of
 * expressions involving vectors and linear operators. As an example
 * consider the addition of multiple vectors
 * @code
 *   dealii::Vector<double> a, b, c, d;
 *   // ..
 *   dealii::Vector<double> result = a + b - c + d;
 * @endcode
 * or the computation of a residual $b-Ax$:
 * @code
 *   dealii::SparseMatrix<double> A;
 *   dealii::Vector<double> b, x;
 *   // ..
 *   const auto op_a = linear_operator(A);
 *
 *   const auto residual =  b - op_a * x;
 * @endcode
 * The expression <code>residual</code> is of type
 * <code>Computation<dealii::Vector<double>></code>. It stores references
 * to <code>A</code>, <code>b</code> and <code>x</code> and defers the
 * actual computation until <code>result</code>, or <code>result_add</code>
 * are explicitly invoked,
 * @code
 *   dealii::Vector<double> y;
 *   residual.reinit_vector(y);
 *   residual.result(y);
 *   residual.result_add(y);
 * @endcode
 * or until the @p Computation object is implictly converted:
 * @code
 *   dealii::Vector<double> y;
 *   y = residual;
 *   y += residual;
 *   y -= residual;
 * @endcode
 *
 * @author Matthias Maier, 2015
 *
 * @ingroup LAOperators
 */
template <typename Range> class Computation
{
public:

  /**
   * Create an empty Computation object. All <code>std::function</code>
   * member objects are initialized with default variants that throw an
   * exception upon invocation.
   */
  Computation ()
  {
    result = [](Range &)
    {
      Assert(false,
             ExcMessage("Uninitialized Computation<Range>::result called"));
    };

    result_add = [](Range &)
    {
      Assert(false,
             ExcMessage("Uninitialized Computation<Range>::result_add called"));
    };

    reinit_vector = [](Range &, bool)
    {
      Assert(false,
             ExcMessage("Uninitialized Computation<Range>::reinit_vector "
                        "method called"));
    };
  }

  /**
   * Default copy constructor.
   */
  Computation (const Computation<Range> &) = default;

  /**
   * Constructor that creates a Computation object from a reference vector
   * @p u. The Computation returns @p u.
   *
   * The Computation object that is created stores a reference to @p u.
   * Thus, the vector must remain a valid reference for the whole lifetime
   * of the Computation object. All changes made on @p u after the creation
   * of the Computation object are reflected by the operator object.
   */
  Computation (const Range &u)
  {
    *this = u;
  }

  /**
   * Default copy assignment operator.
   */
  Computation<Range> &operator=(const Computation<Range> &) = default;

  /**
   * Copy assignment operator that creates a Computation object from a
   * reference vector @p u. The Computation returns @p u.
   *
   * The Computation object that is created stores a reference to @p u.
   * Thus, the vector must remain a valid reference for the whole lifetime
   * of the Computation object. All changes made on @p u after the creation
   * of the Computation object are reflected by the operator object.
   */
  Computation<Range> &operator=(const Range &u)
  {
    result = [&u](Range &v)
    {
      v = u;
    };

    result_add = [&u](Range &v)
    {
      v += u;
    };

    reinit_vector = [&u](Range &v, bool fast)
    {
      v.reinit(u, fast);
    };

    return *this;
  }

  /**
   * Convert a Computation to its result.
   *
   * This conversion operator creates a vector of the Range space and calls
   * <code>result()</code> on it.
   */
  operator Range() const
  {
    Range result_vector;

    reinit_vector(result_vector, /*bool fast=*/ true);
    result(result_vector);

    return result_vector;
  }

  /**
   * @name In-place vector space operations
   */
  //@{

  /**
   * Addition with a Computation @p second_comp with the same @p Range.
   */
  Computation<Range> &operator+=(const Computation<Range> &second_comp)
  {
    return *this = *this + second_comp;
  }

  /**
   * Subtraction with a Computation @p second_comp with the same @p Range.
   */
  Computation<Range> &operator-=(const Computation<Range> &second_comp)
  {
    return *this = *this - second_comp;
  }

  /**
   * Add a constant @p offset (of the @p Range space) to the result of a
   * computation.
   */
  Computation<Range> &operator+=(const Range &offset)
  {
    return *this = *this + Computation<Range>(offset);
  }

  /**
   * Subract a constant @p offset (of the @p Range space) from the result
   * of a computation.
   */
  Computation<Range> &operator-=(const Range &offset)
  {
    return *this = *this - Computation<Range>(offset);
  }

  /**
   * Scalar multiplication of the Computation with a @p number.
   */
  Computation<Range> &operator*=(typename Range::value_type number)
  {
    return *this = *this * number;
  }
  //@}

  /**
   * Store the result of the computation in a vector v of the @p Range
   * space.
   */
  std::function<void(Range &v)> result;

  /**
   * Add the result of the computation to a vector v of the @p Range space.
   */
  std::function<void(Range &v)> result_add;

  /**
   * Initializes a vector v of the Range space to be directly usable as the
   * destination parameter in an application of result, or result_add.
   * Similar to the reinit functions of the vector classes, the boolean
   * determines whether a fast initialization is done, i.e., if it is set
   * to false the content of the vector is set to 0.
   */
  std::function<void(Range &v, bool fast)> reinit_vector;
};


/**
 * @name Vector space operations
 */
//@{

/**
 * @relates Computation
 *
 * Addition of two Computation objects @p first_comp and @p second_comp given by
 * vector space addition of the corresponding results.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range>
operator+(const Computation<Range> &first_comp,
          const Computation<Range> &second_comp)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid computation objects by catching first_comp and
  // second_comp by value

  return_comp.result = [first_comp, second_comp](Range &v)
  {
    first_comp.result(v);
    second_comp.result_add(v);
  };

  return_comp.result_add = [first_comp, second_comp](Range &v)
  {
    first_comp.result_add(v);
    second_comp.result_add(v);
  };

  return return_comp;
}

/**
 * @relates Computation
 *
 * Subtraction of two Computation objects @p first_comp and @p second_comp
 * given by vector space addition of the corresponding results.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range>
operator-(const Computation<Range> &first_comp,
          const Computation<Range> &second_comp)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid computation objects by catching first_comp and
  // second_comp by value

  return_comp.result = [first_comp, second_comp](Range &v)
  {
    second_comp.result(v);
    v *= -1.;
    first_comp.result_add(v);
  };

  return_comp.result_add = [first_comp, second_comp](Range &v)
  {
    first_comp.result_add(v);
    v *= -1.;
    second_comp.result_add(v);
    v *= -1.;
  };

  return return_comp;
}

/**
 * @relates Computation
 *
 * Scalar multiplication of a Computation objects @p comp with a scalar @p number
 * given by a scaling Computation result with @p number.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range>
operator*(const Computation<Range> &comp,
          typename Range::value_type number)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = comp.reinit_vector;

  return_comp.result = [comp, number](Range &v)
  {
    comp.result(v);
    v *= number;
  };

  return_comp.result_add = [comp, number](Range &v)
  {
    v /= number;
    comp.result_add(v);
    v *= number;
  };

  return return_comp;
}

/**
 * @relates Computation
 *
 * Scalar multiplication of a Computation objects @p comp with a scalar @p number
 * given by a scaling Computation result with @p number.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range>
operator*(typename Range::value_type number,
          const Computation<Range> &comp)
{
  return comp * number;
}

/**
 * @relates Computation
 *
 * Add a constant @p offset (of the @p Range space) to the result of a
 * Computation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range> operator+(const Computation<Range> &comp,
                             const Range &offset)
{
  return comp + Computation<Range>(offset);
}

/**
 * @relates Computation
 *
 * Add a constant @p offset (of the @p Range space) to the result of a
 * Computation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range> operator+(const Range &offset,
                             const Computation<Range> &comp)
{
  return Computation<Range>(offset) + comp;
}

/**
 * @relates Computation
 *
 * Subtract a constant @p offset (of the @p Range space) from the result of a
 * Computation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range> operator-(const Computation<Range> &comp,
                             const Range &offset)
{
  return comp - Computation<Range>(offset);
}


/**
 * @relates Computation
 *
 * Subtract a computational result from a constant @p offset (of the @p
 * Range space). The result is a Computation object that applies this
 * computation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
Computation<Range> operator-(const Range &offset,
                             const Computation<Range> &comp)
{
  return Computation<Range>(offset) - comp;
}

//@}


/**
 * @name Creation of a Computation object
 */
//@{

namespace
{
  // Poor man's trait class that determines whether type T is a vector:
  // FIXME: Implement this as a proprer type trait - similar to
  // isBlockVector

  template <typename T>
  class has_vector_interface
  {
    template <typename C>
    static std::false_type test(...);

    template <typename C>
    static std::true_type test(decltype(&C::operator+=),
                               decltype(&C::operator-=),
                               decltype(&C::l2_norm));

  public:
    // type is std::true_type if Matrix provides vmult_add and Tvmult_add,
    // otherwise it is std::false_type

    typedef decltype(test<T>(0, 0, 0)) type;
  };
}

/**
 * @relates Computation
 *
 * Create a Computation object that stores the addition of two vectors.
 *
 * The Computation object that is created stores a reference to @p u and @p
 * v. Thus, the vectors must remain valid references for the whole lifetime
 * of the Computation object. All changes made on @p u or @p v after the
 * creation of the Computation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */

// FIXME: Only implement specialized variant for actual Vector types
template <typename Range,
          typename = typename std::enable_if<has_vector_interface<Range>::type::value>::type>
Computation<Range> operator+(const Range &u, const Range &v)
{
  Computation<Range> return_comp;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.reinit_vector = [&u](Range &x, bool fast)
  {
    x.reinit(u, fast);
  };

  return_comp.result = [&u, &v](Range &x)
  {
    x = u;
    x += v;
  };

  return_comp.result_add = [&u, &v](Range &x)
  {
    x += u;
    x += v;
  };

  return return_comp;
}


/**
 * @relates Computation
 *
 * Create a Computation object that stores the addition of two vectors.
 *
 * The Computation object that is created stores a reference to @p u and @p
 * v. Thus, the vectors must remain valid references for the whole lifetime
 * of the Computation object. All changes made on @p u or @p v after the
 * creation of the Computation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */

// FIXME: Only implement specialized variant for actual Vector types
template <typename Range,
          typename = typename std::enable_if<has_vector_interface<Range>::type::value>::type>
Computation<Range> operator-(const Range &u, const Range &v)
{
  Computation<Range> return_comp;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.reinit_vector = [&u](Range &x, bool fast)
  {
    x.reinit(u, fast);
  };

  return_comp.result = [&u, &v](Range &x)
  {
    x = u;
    x -= v;
  };

  return_comp.result_add = [&u, &v](Range &x)
  {
    x += u;
    x -= v;
  };

  return return_comp;
}


/**
 * @relates Computation
 *
 * Create a Computation object from a LinearOperator and a reference to a
 * vector @p u of the Domain space. The object stores the computation
 * $\text{op} \,u$ (in matrix notation). <code>return</code>
 * (<code>return_add</code>) are implemented with <code>vmult(__1,u)</code>
 * (<code>vmult_add(__1,u)</code>).
 *
 * The Computation object that is created stores a reference to @p u.
 * Thus, the vector must remain a valid reference for the whole lifetime
 * of the Computation object. All changes made on @p u after the creation
 * of the Computation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
Computation<Range>
operator*(const LinearOperator<Range, Domain> &op,
          const Domain &u)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.result = [op, &u](Range &v)
  {
    op.vmult(v, u);
  };

  return_comp.result_add = [op, &u](Range &v)
  {
    op.vmult_add(v, u);
  };

  return return_comp;
}


/**
 * @relates Computation
 *
 * Create a Computation object from a LinearOperator and a reference to a
 * vector @p u of the Range space. The object stores the computation
 * $\text{op}^T \,u$ (in matrix notation). <code>return</code>
 * (<code>return_add</code>) are implemented with <code>Tvmult(__1,u)</code>
 * (<code>Tvmult_add(__1,u)</code>).
 *
 * The Computation object that is created stores a reference to @p u.
 * Thus, the vector must remain a valid reference for the whole lifetime
 * of the Computation object. All changes made on @p u after the creation
 * of the Computation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
Computation<Domain>
operator*(const Range &u,
          const LinearOperator<Range, Domain> &op)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.result = [op, &u](Domain &v)
  {
    op.Tvmult(v, u);
  };

  return_comp.result_add = [op, &u](Domain &v)
  {
    op.Tvmult_add(v, u);
  };

  return return_comp;
}


/**
 * @relates Computation
 *
 * Composition of a Computation object with a LinearOperator.
 * The object stores the computation $\text{op} \,comp$ (in matrix notation).
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
Computation<Range>
operator*(const LinearOperator<Range, Domain> &op,
          const Computation<Domain> &comp)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.result = [op, comp](Domain &v)
  {
    static GrowingVectorMemory<Range> vector_memory;

    Range *i = vector_memory.alloc();
    op.reinit_domain_vector(*i, /*bool fast =*/ true);

    comp.result(*i);
    op.vmult(v, *i);

    vector_memory.free(i);
  };

  return_comp.result_add = [op, comp](Domain &v)
  {
    static GrowingVectorMemory<Range> vector_memory;

    Range *i = vector_memory.alloc();
    op.reinit_range_vector(*i, /*bool fast =*/ true);

    comp.result(*i);
    op.vmult_add(v, *i);

    vector_memory.free(i);
  };

  return return_comp;
}


/**
 * @relates Computation
 *
 * Composition of a Computation object with a LinearOperator.
 * The object stores the computation $\text{op}^T \,comp$ (in matrix notation).
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain>
Computation<Domain>
operator*(const Computation<Range> &comp,
          const LinearOperator<Range, Domain> &op)
{
  Computation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid computation objects by catching op by value
  // u is catched by reference

  return_comp.result = [op, comp](Domain &v)
  {
    static GrowingVectorMemory<Range> vector_memory;

    Range *i = vector_memory.alloc();
    op.reinit_range_vector(*i, /*bool fast =*/ true);

    comp.result(*i);
    op.Tvmult(v, *i);

    vector_memory.free(i);
  };

  return_comp.result_add = [op, comp](Domain &v)
  {
    static GrowingVectorMemory<Range> vector_memory;

    Range *i = vector_memory.alloc();
    op.reinit_range_vector(*i, /*bool fast =*/ true);

    comp.result(*i);
    op.Tvmult_add(v, *i);

    vector_memory.free(i);
  };

  return return_comp;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CXX11
#endif
