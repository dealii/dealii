// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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

#ifndef dealii_packaged_operation_h
#define dealii_packaged_operation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector_memory.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
template <typename Number>
class Vector;
template <typename Range, typename Domain, typename Payload>
class LinearOperator;
template <typename Range = Vector<double>>
class PackagedOperation;
#endif


/**
 * A class to store a computation.
 *
 * The PackagedOperation class allows lazy evaluation of expressions involving
 * vectors and linear operators. This is done by storing the computational
 * expression and only performing the computation when either the object is
 * implicitly converted to a vector object, or <code>apply</code> (or
 * <code>apply_add</code>) is invoked by hand. This avoids unnecessary
 * temporary storage of intermediate results.
 *
 * The class essentially consists of <code>std::function</code> objects that
 * store the knowledge of how to generate the result of a computation and
 * store it in a vector:
 * @code
 *   std::function<void(Range &)> apply;
 *   std::function<void(Range &)> apply_add;
 * @endcode
 *
 * Similar to the LinearOperator class it also has knowledge about how to
 * initialize a vector of the @p Range space:
 * @code
 *   std::function<void(Range &, bool)> reinit_vector;
 * @endcode
 *
 * As an example consider the addition of multiple vectors
 * @code
 *   Vector<double> a, b, c, d;
 *   // ..
 *   Vector<double> result = a + b - c + d;
 * @endcode
 * or the computation of a residual $b-Ax$:
 * @code
 *   SparseMatrix<double> A;
 *   Vector<double> b, x;
 *   // ..
 *   const auto op_a = linear_operator(A);
 *
 *   auto residual =  b - op_a * x;
 * @endcode
 * The expression <code>residual</code> is of type
 * <code>PackagedOperation<Vector<double>></code>. It stores
 * references to <code>A</code>, <code>b</code> and <code>x</code> and defers
 * the actual computation until <code>apply</code>, or <code>apply_add</code>
 * are explicitly invoked,
 * @code
 *   Vector<double> y;
 *   residual.reinit_vector(y);
 *   residual.apply(y);
 *   residual.apply_add(y);
 * @endcode
 * or until the @p PackagedOperation object is implicitly converted:
 * @code
 *   Vector<double> y;
 *   y = residual;
 *   y += residual;
 *   y -= residual;
 * @endcode
 *
 * @note The step-20 tutorial program has a detailed usage example of the
 * LinearOperator class.
 *
 * @author Matthias Maier, 2015
 *
 * @ingroup LAOperators
 */
template <typename Range>
class PackagedOperation
{
public:
  /**
   * Create an empty PackagedOperation object. All <code>std::function</code>
   * member objects are initialized with default variants that throw an
   * exception upon invocation.
   */
  PackagedOperation()
  {
    apply = [](Range &) {
      Assert(false,
             ExcMessage(
               "Uninitialized PackagedOperation<Range>::apply called"));
    };

    apply_add = [](Range &) {
      Assert(false,
             ExcMessage(
               "Uninitialized PackagedOperation<Range>::apply_add called"));
    };

    reinit_vector = [](Range &, bool) {
      Assert(false,
             ExcMessage("Uninitialized PackagedOperation<Range>::reinit_vector "
                        "method called"));
    };
  }

  /**
   * Default copy constructor.
   */
  PackagedOperation(const PackagedOperation<Range> &) = default;

  /**
   * Constructor that creates a PackagedOperation object from a reference
   * vector @p u. The PackagedOperation returns @p u.
   *
   * The PackagedOperation object that is created stores a reference to @p u.
   * Thus, the vector must remain a valid reference for the whole lifetime of
   * the PackagedOperation object. All changes made on @p u after the creation
   * of the PackagedOperation object are reflected by the operator object.
   */
  PackagedOperation(const Range &u)
  {
    *this = u;
  }

  /**
   * Default copy assignment operator.
   */
  PackagedOperation<Range> &
  operator=(const PackagedOperation<Range> &) = default;

  /**
   * Copy assignment operator that creates a PackagedOperation object from a
   * reference vector @p u. The PackagedOperation returns @p u.
   *
   * The PackagedOperation object that is created stores a reference to @p u.
   * Thus, the vector must remain a valid reference for the whole lifetime of
   * the PackagedOperation object. All changes made on @p u after the creation
   * of the PackagedOperation object are reflected by the operator object.
   */
  PackagedOperation<Range> &
  operator=(const Range &u)
  {
    apply = [&u](Range &v) { v = u; };

    apply_add = [&u](Range &v) { v += u; };

    reinit_vector = [&u](Range &v, bool omit_zeroing_entries) {
      v.reinit(u, omit_zeroing_entries);
    };

    return *this;
  }

  /**
   * Convert a PackagedOperation to its result.
   *
   * This conversion operator creates a vector of the Range space and calls
   * <code>apply()</code> on it.
   */
  operator Range() const
  {
    Range result_vector;

    reinit_vector(result_vector, /*bool omit_zeroing_entries=*/true);
    apply(result_vector);

    return result_vector;
  }

  /**
   * @name In-place vector space operations
   */
  //@{

  /**
   * Addition with a PackagedOperation @p second_comp with the same @p Range.
   */
  PackagedOperation<Range> &
  operator+=(const PackagedOperation<Range> &second_comp)
  {
    *this = *this + second_comp;
    return *this;
  }

  /**
   * Subtraction with a PackagedOperation @p second_comp with the same @p
   * Range.
   */
  PackagedOperation<Range> &
  operator-=(const PackagedOperation<Range> &second_comp)
  {
    *this = *this - second_comp;
    return *this;
  }

  /**
   * Add a constant @p offset (of the @p Range space) to the result of a
   * PackagedOperation.
   */
  PackagedOperation<Range> &
  operator+=(const Range &offset)
  {
    *this = *this + PackagedOperation<Range>(offset);
    return *this;
  }

  /**
   * Subtract a constant @p offset (of the @p Range space) from the result of
   * a PackagedOperation.
   */
  PackagedOperation<Range> &
  operator-=(const Range &offset)
  {
    *this = *this - PackagedOperation<Range>(offset);
    return *this;
  }

  /**
   * Scalar multiplication of the PackagedOperation with a @p number.
   */
  PackagedOperation<Range> &
  operator*=(typename Range::value_type number)
  {
    *this = *this * number;
    return *this;
  }
  //@}

  /**
   * Store the result of the PackagedOperation in a vector v of the @p Range
   * space.
   */
  std::function<void(Range &v)> apply;

  /**
   * Add the result of the PackagedOperation to a vector v of the @p Range
   * space.
   */
  std::function<void(Range &v)> apply_add;

  /**
   * Initializes a vector v of the Range space to be directly usable as the
   * destination parameter in an application of apply, or apply_add. Similar
   * to the reinit functions of the vector classes, the boolean determines
   * whether a fast initialization is done, i.e., if it is set to false the
   * content of the vector is set to 0.
   */
  std::function<void(Range &v, bool omit_zeroing_entries)> reinit_vector;
};


/**
 * @name Vector space operations
 */
//@{

/**
 * @relatesalso PackagedOperation
 *
 * Addition of two PackagedOperation objects @p first_comp and @p second_comp
 * given by vector space addition of the corresponding results.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator+(const PackagedOperation<Range> &first_comp,
          const PackagedOperation<Range> &second_comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid PackagedOperation objects by catching first_comp and
  // second_comp by value

  return_comp.apply = [first_comp, second_comp](Range &v) {
    first_comp.apply(v);
    second_comp.apply_add(v);
  };

  return_comp.apply_add = [first_comp, second_comp](Range &v) {
    first_comp.apply_add(v);
    second_comp.apply_add(v);
  };

  return return_comp;
}

/**
 * @relatesalso PackagedOperation
 *
 * Subtraction of two PackagedOperation objects @p first_comp and @p
 * second_comp given by vector space addition of the corresponding results.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator-(const PackagedOperation<Range> &first_comp,
          const PackagedOperation<Range> &second_comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid PackagedOperation objects by catching first_comp and
  // second_comp by value

  return_comp.apply = [first_comp, second_comp](Range &v) {
    second_comp.apply(v);
    v *= -1.;
    first_comp.apply_add(v);
  };

  return_comp.apply_add = [first_comp, second_comp](Range &v) {
    first_comp.apply_add(v);
    v *= -1.;
    second_comp.apply_add(v);
    v *= -1.;
  };

  return return_comp;
}

/**
 * @relatesalso PackagedOperation
 *
 * Scalar multiplication of a PackagedOperation objects @p comp with a scalar
 * @p number given by a scaling PackagedOperation result with @p number.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range> operator*(const PackagedOperation<Range> &comp,
                                   typename Range::value_type      number)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = comp.reinit_vector;

  // the trivial case: number is zero
  if (number == 0.)
    {
      return_comp.apply = [](Range &v) { v = 0.; };

      return_comp.apply_add = [](Range &) {};
    }
  else
    {
      return_comp.apply = [comp, number](Range &v) {
        comp.apply(v);
        v *= number;
      };

      return_comp.apply_add = [comp, number](Range &v) {
        v /= number;
        comp.apply_add(v);
        v *= number;
      };
    }

  return return_comp;
}

/**
 * @relatesalso PackagedOperation
 *
 * Scalar multiplication of a PackagedOperation objects @p comp with a scalar
 * @p number given by a scaling PackagedOperation result with @p number.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range> operator*(typename Range::value_type      number,
                                   const PackagedOperation<Range> &comp)
{
  return comp * number;
}

/**
 * @relatesalso PackagedOperation
 *
 * Add a constant @p offset (of the @p Range space) to the result of a
 * PackagedOperation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator+(const PackagedOperation<Range> &comp, const Range &offset)
{
  return comp + PackagedOperation<Range>(offset);
}

/**
 * @relatesalso PackagedOperation
 *
 * Add a constant @p offset (of the @p Range space) to the result of a
 * PackagedOperation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator+(const Range &offset, const PackagedOperation<Range> &comp)
{
  return PackagedOperation<Range>(offset) + comp;
}

/**
 * @relatesalso PackagedOperation
 *
 * Subtract a constant @p offset (of the @p Range space) from the result of a
 * PackagedOperation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator-(const PackagedOperation<Range> &comp, const Range &offset)
{
  return comp - PackagedOperation<Range>(offset);
}


/**
 * @relatesalso PackagedOperation
 *
 * Subtract a computational result from a constant @p offset (of the @p Range
 * space). The result is a PackagedOperation object that applies this
 * computation.
 *
 * @ingroup LAOperators
 */
template <typename Range>
PackagedOperation<Range>
operator-(const Range &offset, const PackagedOperation<Range> &comp)
{
  return PackagedOperation<Range>(offset) - comp;
}

//@}


/**
 * @name Creation of a PackagedOperation object
 */
//@{

namespace internal
{
  namespace PackagedOperationImplementation
  {
    // Poor man's trait class that determines whether type T is a vector:
    // FIXME: Implement this as a proper type trait - similar to
    // isBlockVector

    template <typename T>
    class has_vector_interface
    {
      template <typename C>
      static std::false_type
      test(...);

      template <typename C>
      static std::true_type
      test(decltype(&C::operator+=),
           decltype(&C::operator-=),
           decltype(&C::l2_norm));

    public:
      // type is std::true_type if Matrix provides vmult_add and Tvmult_add,
      // otherwise it is std::false_type

      using type = decltype(test<T>(nullptr, nullptr, nullptr));
    }; // namespace
  }    // namespace PackagedOperationImplementation
} // namespace internal


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object that stores the addition of two vectors.
 *
 * The PackagedOperation object that is created stores a reference to @p u and
 * @p v. Thus, the vectors must remain valid references for the whole lifetime
 * of the PackagedOperation object. All changes made on @p u or @p v after the
 * creation of the PackagedOperation object are reflected by the operator
 * object.
 *
 * @ingroup LAOperators
 */

template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range>
operator+(const Range &u, const Range &v)
{
  PackagedOperation<Range> return_comp;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.reinit_vector = [&u](Range &x, bool omit_zeroing_entries) {
    x.reinit(u, omit_zeroing_entries);
  };

  return_comp.apply = [&u, &v](Range &x) {
    x = u;
    x += v;
  };

  return_comp.apply_add = [&u, &v](Range &x) {
    x += u;
    x += v;
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object that stores the subtraction of two
 * vectors.
 *
 * The PackagedOperation object that is created stores a reference to @p u and
 * @p v. Thus, the vectors must remain valid references for the whole lifetime
 * of the PackagedOperation object. All changes made on @p u or @p v after the
 * creation of the PackagedOperation object are reflected by the operator
 * object.
 *
 * @ingroup LAOperators
 */

template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range>
operator-(const Range &u, const Range &v)
{
  PackagedOperation<Range> return_comp;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.reinit_vector = [&u](Range &x, bool omit_zeroing_entries) {
    x.reinit(u, omit_zeroing_entries);
  };

  return_comp.apply = [&u, &v](Range &x) {
    x = u;
    x -= v;
  };

  return_comp.apply_add = [&u, &v](Range &x) {
    x += u;
    x -= v;
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object that stores the scaling of a vector with
 * a @p number.
 *
 * The PackagedOperation object that is created stores a reference to @p u.
 * Thus, the vectors must remain valid references for the whole lifetime of
 * the PackagedOperation object. All changes made on @p u or @p v after the
 * creation of the PackagedOperation object are reflected by the operator
 * object.
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range> operator*(const Range &              u,
                                   typename Range::value_type number)
{
  return PackagedOperation<Range>(u) * number;
}


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object that stores the scaling of a vector with
 * a @p number.
 *
 * The PackagedOperation object that is created stores a reference to @p u.
 * Thus, the vectors must remain valid references for the whole lifetime of
 * the PackagedOperation object. All changes made on @p u or @p v after the
 * creation of the PackagedOperation object are reflected by the operator
 * object.
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range> operator*(typename Range::value_type number,
                                   const Range &              u)
{
  return number * PackagedOperation<Range>(u);
}


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object from a LinearOperator and a reference to
 * a vector @p u of the Domain space. The object stores the PackagedOperation
 * $\text{op} \,u$ (in matrix notation). <code>return</code>
 * (<code>return_add</code>) are implemented with <code>vmult(__1,u)</code>
 * (<code>vmult_add(__1,u)</code>).
 *
 * The PackagedOperation object that is created stores a reference to @p u.
 * Thus, the vector must remain a valid reference for the whole lifetime of
 * the PackagedOperation object. All changes made on @p u after the creation
 * of the PackagedOperation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
operator*(const LinearOperator<Range, Domain, Payload> &op, const Domain &u)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, &u](Range &v) { op.vmult(v, u); };

  return_comp.apply_add = [op, &u](Range &v) { op.vmult_add(v, u); };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation
 *
 * Create a PackagedOperation object from a LinearOperator and a reference to
 * a vector @p u of the Range space. The object stores the PackagedOperation
 * $\text{op}^T \,u$ (in matrix notation). <code>return</code>
 * (<code>return_add</code>) are implemented with <code>Tvmult(__1,u)</code>
 * (<code>Tvmult_add(__1,u)</code>).
 *
 * The PackagedOperation object that is created stores a reference to @p u.
 * Thus, the vector must remain a valid reference for the whole lifetime of
 * the PackagedOperation object. All changes made on @p u after the creation
 * of the PackagedOperation object are reflected by the operator object.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Domain>
operator*(const Range &u, const LinearOperator<Range, Domain, Payload> &op)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, &u](Domain &v) { op.Tvmult(v, u); };

  return_comp.apply_add = [op, &u](Domain &v) { op.Tvmult_add(v, u); };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation
 *
 * Composition of a PackagedOperation object with a LinearOperator. The object
 * stores the computation $\text{op} \,comp$ (in matrix notation).
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
operator*(const LinearOperator<Range, Domain, Payload> &op,
          const PackagedOperation<Domain> &             comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_domain_vector(*i, /*bool omit_zeroing_entries =*/true);

    comp.apply(*i);
    op.vmult(v, *i);
  };

  return_comp.apply_add = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i, /*bool omit_zeroing_entries =*/true);

    comp.apply(*i);
    op.vmult_add(v, *i);
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation
 *
 * Composition of a PackagedOperation object with a LinearOperator. The object
 * stores the computation $\text{op}^T \,comp$ (in matrix notation).
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Domain>
operator*(const PackagedOperation<Range> &              comp,
          const LinearOperator<Range, Domain, Payload> &op)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i, /*bool omit_zeroing_entries =*/true);

    comp.apply(*i);
    op.Tvmult(v, *i);
  };

  return_comp.apply_add = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i, /*bool omit_zeroing_entries =*/true);

    comp.apply(*i);
    op.Tvmult_add(v, *i);
  };

  return return_comp;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif
