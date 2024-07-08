// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_linear_operator_h
#define dealii_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector_memory.h>

#include <array>
#include <functional>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
namespace internal
{
  namespace LinearOperatorImplementation
  {
    class EmptyPayload;
  }
} // namespace internal

template <typename Number>
class Vector;

class PreconditionIdentity;

template <typename Range  = Vector<double>,
          typename Domain = Range,
          typename Payload =
            internal::LinearOperatorImplementation::EmptyPayload>
class LinearOperator;
#endif

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &, const Matrix &);

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const Matrix &);

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Domain, Payload>
null_operator(const LinearOperator<Range, Domain, Payload> &);

template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
identity_operator(const LinearOperator<Range, Domain, Payload> &);


/**
 * A class to store the abstract concept of a linear operator.
 *
 * The class essentially consists of <code>std::function</code> objects that
 * store the knowledge of how to apply the linear operator by implementing the
 * abstract @p Matrix interface:
 * @code
 *   std::function<void(Range &, const Domain &)> vmult;
 *   std::function<void(Range &, const Domain &)> vmult_add;
 *   std::function<void(Domain &, const Range &)> Tvmult;
 *   std::function<void(Domain &, const Range &)> Tvmult_add;
 * @endcode
 *
 * But, in contrast to a usual matrix object, the domain and range of the
 * linear operator are also bound to the LinearOperator class on the type
 * level. Because of this, `LinearOperator<Range, Domain>` has two
 * additional function objects
 * @code
 *   std::function<void(Range &, bool)> reinit_range_vector;
 *   std::function<void(Domain &, bool)> reinit_domain_vector;
 * @endcode
 * that store the knowledge how to initialize (resize + internal data
 * structures) an arbitrary vector of the @p Range and @p Domain space.
 *
 * The primary purpose of this class is to provide syntactic sugar for complex
 * matrix-vector operations and free the user from having to create, set up
 * and handle intermediate storage locations by hand.
 *
 * As an example consider the operation $(A+k\,B)\,C$, where $A$, $B$ and $C$
 * denote (possible different) matrices. In order to construct a
 * LinearOperator <code>op</code> that stores the knowledge of this operation,
 * one can write:
 *
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 *
 * dealii::SparseMatrix<double> A, B, C;
 * const double k = ...;
 *
 * // Setup and assembly of matrices
 *
 * const auto op_a = linear_operator(A);
 * const auto op_b = linear_operator(B);
 * const auto op_c = linear_operator(C);
 *
 * const auto op = (op_a + k * op_b) * op_c;
 * @endcode
 *
 * @note This class makes heavy use of <code>std::function</code> objects and
 * lambda functions. This flexibility comes with a run-time penalty. Only use
 * this object to encapsulate matrix object of medium to large size (as a rule
 * of thumb, sparse matrices with a size $1000\times1000$, or larger).
 *
 * @note In order to use Trilinos or PETSc sparse matrices and preconditioners
 * in conjunction with the LinearOperator class, it is necessary to extend the
 * functionality of the LinearOperator class by means of an additional Payload.
 *
 * For example: LinearOperator instances representing matrix inverses usually
 * require calling some linear solver. These solvers may not have interfaces
 * to the LinearOperator (which, for example, may represent a composite
 * operation). The
 * TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
 * therefore provides an interface extension to the LinearOperator so that it
 * can be passed to the solver and used by the solver as if it were a Trilinos
 * operator. This implies that all of the necessary functionality of the
 * specific Trilinos operator has been overloaded within the Payload class.
 * This includes operator-vector multiplication and inverse operator-vector
 * multiplication, where the operator can be either a
 * TrilinosWrappers::SparseMatrix or a TrilinosWrappers::PreconditionBase
 * and the vector is a native Trilinos vector.
 *
 * Another case where payloads provide a crucial supplement to the
 * LinearOperator class are when composite operations are constructed (via
 * operator overloading). In this instance, it is again necessary to provide
 * an interface that produces the result of this composite operation that is
 * compatible with Trilinos operator used by Trilinos solvers.
 *
 * @note Many use cases of LinearOperator lead to intermediate expressions
 * requiring a PackagedOperation. In order to include all necessary header
 * files in one go consider using
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 * @endcode
 *
 * In order to use the full LinearOperator and PackagedOperation
 *
 * @note To ensure that the correct payload is provided, wrapper functions
 * for linear operators have been provided within the respective
 * TrilinosWrappers (and, in the future, PETScWrappers) namespaces.
 *
 * <h3> Examples of use </h3>
 * The step-20 tutorial program has a detailed usage example of the
 * LinearOperator class.
 *
 * <h3> Instrumenting operations </h3>
 * It is sometimes useful to know when functions are called, or to inject
 * additional operations. In such cases, what one wants is to replace, for
 * example, the `vmult` object of this class with one that does the additional
 * operations and then calls what was originally supposed to happen. This
 * can be done with commands such as the following:
 * @code
 *    auto A_inv  = inverse_operator(A, solver_A, preconditioner_A);
 *    A_inv.vmult = [base_vmult = A_inv.vmult](Vector<double>       &dst,
 *                                             const Vector<double> &src) {
 *      std::cout << "Calling A_inv.vmult()" << std::endl;
 *      base_vmult(dst, src);
 *    };
 * @endcode
 * Here, we replace `A_inv.vmult` with a lambda function that first captures
 * the previous value of `A_inv.vmult` and stores it in the `base_vmult`
 * object. The newly installed `A_inv.vmult` function then first outputs some
 * status information, and then calls the original functionality.
 *
 * This approach works for all of the other function objects mentioned above
 * as well.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
class LinearOperator : public Payload
{
public:
  /**
   * Create an empty LinearOperator object.
   * When a payload is passed to this constructor, the resulting operator is
   * constructed with a functional payload.
   * In either case, this constructor yields an object that can not actually
   * be used for any linear operator operations, and will throw an exception
   * upon invocation.
   */
  LinearOperator(const Payload &payload = Payload())
    : Payload(payload)
    , is_null_operator(false)
  {
    vmult = [](Range &, const Domain &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::vmult called"));
    };

    vmult_add = [](Range &, const Domain &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::vmult_add called"));
    };

    Tvmult = [](Domain &, const Range &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::Tvmult called"));
    };

    Tvmult_add = [](Domain &, const Range &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::Tvmult_add called"));
    };

    reinit_range_vector = [](Range &, bool) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::reinit_range_vector method called"));
    };

    reinit_domain_vector = [](Domain &, bool) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::reinit_domain_vector method called"));
    };
  }

  /**
   * Default copy constructor.
   */
  LinearOperator(const LinearOperator<Range, Domain, Payload> &) = default;

  /**
   * Templated copy constructor that creates a LinearOperator object from an
   * object @p op for which the conversion function
   * <code>linear_operator</code> is defined.
   */
  template <typename Op,
            typename = std::enable_if_t<
              !std::is_base_of_v<LinearOperator<Range, Domain, Payload>, Op>>>
  LinearOperator(const Op &op)
  {
    *this = linear_operator<Range, Domain, Payload, Op>(op);
  }

  /**
   * Default copy assignment operator.
   */
  LinearOperator<Range, Domain, Payload> &
  operator=(const LinearOperator<Range, Domain, Payload> &) = default;

  /**
   * Templated copy assignment operator for an object @p op for which the
   * conversion function <code>linear_operator</code> is defined.
   */
  template <typename Op,
            typename = std::enable_if_t<
              !std::is_base_of_v<LinearOperator<Range, Domain, Payload>, Op>>>
  LinearOperator<Range, Domain, Payload> &
  operator=(const Op &op)
  {
    *this = linear_operator<Range, Domain, Payload, Op>(op);
    return *this;
  }

  /**
   * Application of the LinearOperator object to a vector u of the @p Domain
   * space giving a vector v of the @p Range space.
   */
  std::function<void(Range &v, const Domain &u)> vmult;

  /**
   * Application of the LinearOperator object to a vector u of the @p Domain
   * space. The result is added to the vector v.
   */
  std::function<void(Range &v, const Domain &u)> vmult_add;

  /**
   * Application of the transpose LinearOperator object to a vector u of the
   * @p Range space giving a vector v of the @p Domain space.
   */
  std::function<void(Domain &v, const Range &u)> Tvmult;

  /**
   * Application of the transpose LinearOperator object to a vector @p u of
   * the @p Range space.The result is added to the vector @p v.
   */
  std::function<void(Domain &v, const Range &u)> Tvmult_add;

  /**
   * Initializes a vector v of the Range space to be directly usable as the
   * destination parameter in an application of vmult. Similar to the reinit
   * functions of the vector classes, the boolean determines whether a fast
   * initialization is done, i.e., if it is set to false the content of the
   * vector is set to 0.
   */
  std::function<void(Range &v, bool omit_zeroing_entries)> reinit_range_vector;

  /**
   * Initializes a vector of the Domain space to be directly usable as the
   * source parameter in an application of vmult. Similar to the reinit
   * functions of the vector classes, the boolean determines whether a fast
   * initialization is done, i.e., if it is set to false the content of the
   * vector is set to 0.
   */
  std::function<void(Domain &v, bool omit_zeroing_entries)>
    reinit_domain_vector;

  /**
   * @name In-place vector space operations
   */
  /** @{ */

  /**
   * Addition with a LinearOperator @p second_op with the same @p Domain and
   * @p Range.
   */
  LinearOperator<Range, Domain, Payload> &
  operator+=(const LinearOperator<Range, Domain, Payload> &second_op)
  {
    *this = *this + second_op;
    return *this;
  }

  /**
   * Subtraction with a LinearOperator @p second_op with the same @p Domain
   * and @p Range.
   */
  LinearOperator<Range, Domain, Payload> &
  operator-=(const LinearOperator<Range, Domain, Payload> &second_op)
  {
    *this = *this - second_op;
    return *this;
  }

  /**
   * Composition of the LinearOperator with an endomorphism @p second_op of
   * the @p Domain space.
   */
  LinearOperator<Range, Domain, Payload> &
  operator*=(const LinearOperator<Domain, Domain, Payload> &second_op)
  {
    *this = *this * second_op;
    return *this;
  }

  /**
   * Scalar multiplication of the LinearOperator with @p number from the
   * right.
   */
  LinearOperator<Range, Domain, Payload>
  operator*=(typename Domain::value_type number)
  {
    *this = *this * number;
    return *this;
  }

  /**
   * This bool is used to determine whether a linear operator is a null
   * operator. In this case the class is able to optimize some operations like
   * multiplication or addition.
   */
  bool is_null_operator;

  /** @} */
};


/**
 * @name Vector space operations
 */
/** @{ */

/**
 * @relatesalso LinearOperator
 *
 * Addition of two linear operators @p first_op and @p second_op given by
 * $(\mathrm{first\_op}+\mathrm{second\_op})x \dealcoloneq \mathrm{first\_op}(x)
 * + \mathrm{second\_op}(x)$
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator+(const LinearOperator<Range, Domain, Payload> &first_op,
          const LinearOperator<Range, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator)
    {
      return second_op;
    }
  else if (second_op.is_null_operator)
    {
      return first_op;
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op{
        static_cast<const Payload &>(first_op) +
        static_cast<const Payload &>(second_op)};

      return_op.reinit_range_vector  = first_op.reinit_range_vector;
      return_op.reinit_domain_vector = first_op.reinit_domain_vector;

      // ensure to have valid computation objects by catching first_op and
      // second_op by value

      return_op.vmult = [first_op, second_op](Range &v, const Domain &u) {
        first_op.vmult(v, u);
        second_op.vmult_add(v, u);
      };

      return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u) {
        first_op.vmult_add(v, u);
        second_op.vmult_add(v, u);
      };

      return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u) {
        second_op.Tvmult(v, u);
        first_op.Tvmult_add(v, u);
      };

      return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u) {
        second_op.Tvmult_add(v, u);
        first_op.Tvmult_add(v, u);
      };

      return return_op;
    }
}


/**
 * @relatesalso LinearOperator
 *
 * Subtraction of two linear operators @p first_op and @p second_op given by
 * $(\mathrm{first\_op}-\mathrm{second\_op})x \dealcoloneq \mathrm{first\_op}(x)
 * - \mathrm{second\_op}(x)$
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator-(const LinearOperator<Range, Domain, Payload> &first_op,
          const LinearOperator<Range, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator)
    {
      return -1. * second_op;
    }
  else if (second_op.is_null_operator)
    {
      return first_op;
    }
  else
    {
      // implement with addition and scalar multiplication
      return first_op + (-1. * second_op);
    }
}


/**
 * @relatesalso LinearOperator
 *
 * Scalar multiplication of a ScalarOperator object @p op with @p number from
 * the left.
 *
 * The @p Domain and @p Range types must implement the following
 * <code>operator*=</code> member functions accepting the appropriate scalar
 * Number type for rescaling:
 *
 * @code
 * Domain & operator *=(Domain::value_type);
 * Range & operator *=(Range::value_type);
 * @endcode
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(typename Range::value_type                    number,
          const LinearOperator<Range, Domain, Payload> &op)
{
  static_assert(
    std::is_convertible_v<typename Range::value_type,
                          typename Domain::value_type>,
    "Range and Domain must have implicitly convertible 'value_type's");

  if (op.is_null_operator)
    {
      return op;
    }
  else if (number == 0.)
    {
      return null_operator(op);
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op = op;

      // ensure to have valid computation objects by catching number and op by
      // value

      return_op.vmult = [number, op](Range &v, const Domain &u) {
        op.vmult(v, u);
        v *= number;
      };

      return_op.vmult_add = [number, op](Range &v, const Domain &u) {
        v /= number;
        op.vmult_add(v, u);
        v *= number;
      };

      return_op.Tvmult = [number, op](Domain &v, const Range &u) {
        op.Tvmult(v, u);
        v *= number;
      };

      return_op.Tvmult_add = [number, op](Domain &v, const Range &u) {
        v /= number;
        op.Tvmult_add(v, u);
        v *= number;
      };

      return return_op;
    }
}


/**
 * @relatesalso LinearOperator
 *
 * Scalar multiplication of a ScalarOperator object from the right.
 *
 * The @p Domain and @p Range types must implement the following
 * <code>operator*=</code> member functions for rescaling:
 *
 * @code
 * Domain & operator *=(Domain::value_type);
 * Range & operator *=(Range::value_type);
 * @endcode
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(const LinearOperator<Range, Domain, Payload> &op,
          typename Domain::value_type                   number)
{
  static_assert(
    std::is_convertible_v<typename Range::value_type,
                          typename Domain::value_type>,
    "Range and Domain must have implicitly convertible 'value_type's");

  return number * op;
}

/** @} */


/**
 * @name Composition and manipulation of a LinearOperator
 */
/** @{ */

/**
 * @relatesalso LinearOperator
 *
 * Composition of two linear operators @p first_op and @p second_op given by
 * $(\mathrm{first\_op}*\mathrm{second\_op})x \dealcoloneq
 * \mathrm{first\_op}(\mathrm{second\_op}(x))$
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename Intermediate,
          typename Domain,
          typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(const LinearOperator<Range, Intermediate, Payload>  &first_op,
          const LinearOperator<Intermediate, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator || second_op.is_null_operator)
    {
      LinearOperator<Range, Domain, Payload> return_op;
      return_op.reinit_domain_vector = second_op.reinit_domain_vector;
      return_op.reinit_range_vector  = first_op.reinit_range_vector;
      return null_operator(return_op);
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op{
        static_cast<const Payload &>(first_op) *
        static_cast<const Payload &>(second_op)};

      return_op.reinit_domain_vector = second_op.reinit_domain_vector;
      return_op.reinit_range_vector  = first_op.reinit_range_vector;

      // ensure to have valid computation objects by catching first_op and
      // second_op by value

      return_op.vmult = [first_op, second_op](Range &v, const Domain &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        second_op.reinit_range_vector(*i, /*bool omit_zeroing_entries =*/true);
        second_op.vmult(*i, u);
        first_op.vmult(v, *i);
      };

      return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        second_op.reinit_range_vector(*i, /*bool omit_zeroing_entries =*/true);
        second_op.vmult(*i, u);
        first_op.vmult_add(v, *i);
      };

      return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        first_op.reinit_domain_vector(*i, /*bool omit_zeroing_entries =*/true);
        first_op.Tvmult(*i, u);
        second_op.Tvmult(v, *i);
      };

      return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        first_op.reinit_domain_vector(*i, /*bool omit_zeroing_entries =*/true);
        first_op.Tvmult(*i, u);
        second_op.Tvmult_add(v, *i);
      };

      return return_op;
    }
}


/**
 * @relatesalso LinearOperator
 *
 * Return the transpose linear operations of @p op.
 *
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Domain, Range, Payload>
transpose_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  LinearOperator<Domain, Range, Payload> return_op{op.transpose_payload()};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult      = op.Tvmult;
  return_op.vmult_add  = op.Tvmult_add;
  return_op.Tvmult     = op.vmult;
  return_op.Tvmult_add = op.vmult_add;

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Return an object representing the inverse of the LinearOperator @p op.
 *
 * The function takes references @p solver and @p preconditioner to an
 * iterative solver and a preconditioner that are used in the
 * <code>vmult</code> and <code>Tvmult</code> implementations of the
 * LinearOperator object.
 *
 * The LinearOperator object that is created stores a reference to @p solver
 * and @p preconditioner. Thus, both objects must remain a valid reference for
 * the whole lifetime of the LinearOperator object. Internal data structures
 * of the @p solver object will be modified upon invocation of
 * <code>vmult</code> or <code>Tvmult</code>.
 *
 *
 * @ingroup LAOperators
 */
template <typename Payload,
          typename Solver,
          typename Preconditioner,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver                                       &solver,
                 const Preconditioner                         &preconditioner)
{
  LinearOperator<Domain, Range, Payload> return_op{
    op.inverse_payload(solver, preconditioner)};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = [op, &solver, &preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v, /*bool omit_zeroing_entries =*/false);
    solver.solve(op, v, u, preconditioner);
  };

  return_op.vmult_add = [op, &solver, &preconditioner](Range        &v,
                                                       const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2, /*bool omit_zeroing_entries =*/false);
    solver.solve(op, *v2, u, preconditioner);
    v += *v2;
  };

  return_op.Tvmult = [op, &solver, &preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v, /*bool omit_zeroing_entries =*/false);
    solver.solve(transpose_operator(op), v, u, preconditioner);
  };

  return_op.Tvmult_add = [op, &solver, &preconditioner](Range        &v,
                                                        const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2, /*bool omit_zeroing_entries =*/false);
    solver.solve(transpose_operator(op), *v2, u, preconditioner);
    v += *v2;
  };

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Variant of above function that takes a LinearOperator @p preconditioner
 * as preconditioner argument.
 *
 * @ingroup LAOperators
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver                                       &solver,
                 const LinearOperator<Range, Domain, Payload> &preconditioner)
{
  LinearOperator<Domain, Range, Payload> return_op{
    op.inverse_payload(solver, preconditioner)};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = [op, &solver, preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v, /*bool omit_zeroing_entries =*/false);
    solver.solve(op, v, u, preconditioner);
  };

  return_op.vmult_add = [op, &solver, preconditioner](Range        &v,
                                                      const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2, /*bool omit_zeroing_entries =*/false);
    solver.solve(op, *v2, u, preconditioner);
    v += *v2;
  };

  return_op.Tvmult = [op, &solver, preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v, /*bool omit_zeroing_entries =*/false);
    solver.solve(transpose_operator(op), v, u, preconditioner);
  };

  return_op.Tvmult_add = [op, &solver, preconditioner](Range        &v,
                                                       const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2, /*bool omit_zeroing_entries =*/false);
    solver.solve(transpose_operator(op), *v2, u, preconditioner);
    v += *v2;
  };

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Variant of above function without a preconditioner argument. In this
 * case the identity_operator() of the @p op argument is used as a
 * preconditioner. This is equivalent to using PreconditionIdentity.
 *
 * @ingroup LAOperators
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver                                       &solver)
{
  return inverse_operator(op, solver, identity_operator(op));
}


/**
 * @relatesalso LinearOperator
 *
 * Special overload of above function that takes a PreconditionIdentity
 * argument.
 *
 * @ingroup LAOperators
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver                                       &solver,
                 const PreconditionIdentity &)
{
  return inverse_operator(op, solver);
}

/** @} */


/**
 * @name Creation of a LinearOperator
 */
/** @{ */

/**
 * @relatesalso LinearOperator
 *
 * Return a LinearOperator that is the identity of the vector space @p Range.
 *
 * The function takes an <code>std::function</code> object @p reinit_vector as
 * an argument to initialize the <code>reinit_range_vector</code> and
 * <code>reinit_domain_vector</code> objects of the LinearOperator object.
 *
 * @ingroup LAOperators
 */
template <
  typename Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Range, Payload>
identity_operator(const std::function<void(Range &, bool)> &reinit_vector)
{
  LinearOperator<Range, Range, Payload> return_op{Payload()};

  return_op.reinit_range_vector  = reinit_vector;
  return_op.reinit_domain_vector = reinit_vector;

  return_op.vmult = [](Range &v, const Range &u) { v = u; };

  return_op.vmult_add = [](Range &v, const Range &u) { v += u; };

  return_op.Tvmult = [](Range &v, const Range &u) { v = u; };

  return_op.Tvmult_add = [](Range &v, const Range &u) { v += u; };

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Return a LinearOperator that is the identity of the vector space @p Range.
 *
 * The function takes a LinearOperator @p op and uses its range initializer
 * to create an identity operator. In contrast to the function above, this
 * function also ensures that the underlying Payload matches that of the
 * input @p op.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
identity_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  auto return_op = identity_operator<Range, Payload>(op.reinit_range_vector);
  static_cast<Payload &>(return_op) = op.identity_payload();

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Return a nulled variant of the LinearOperator @p op, i.e. with optimized
 * LinearOperator::vmult, LinearOperator::vmult_add, etc. functions and with
 * LinearOperator::is_null_operator set to true.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
null_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  LinearOperator<Range, Domain, Payload> return_op{op.null_payload()};

  return_op.is_null_operator = true;

  return_op.reinit_range_vector  = op.reinit_range_vector;
  return_op.reinit_domain_vector = op.reinit_domain_vector;

  return_op.vmult = [](Range &v, const Domain &) { v = 0.; };

  return_op.vmult_add = [](Range &, const Domain &) {};

  return_op.Tvmult = [](Domain &v, const Range &) { v = 0.; };

  return_op.Tvmult_add = [](Domain &, const Range &) {};

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Return a LinearOperator that acts as a mean value filter. The vmult()
 * functions of this matrix subtract the mean values of the vector.
 *
 * The function takes an <code>std::function</code> object @p reinit_vector as
 * an argument to initialize the <code>reinit_range_vector</code> and
 * <code>reinit_domain_vector</code> objects of the LinearOperator object.
 *
 * @ingroup LAOperators
 */
template <
  typename Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Range, Payload>
mean_value_filter(const std::function<void(Range &, bool)> &reinit_vector)
{
  LinearOperator<Range, Range, Payload> return_op{Payload()};

  return_op.reinit_range_vector  = reinit_vector;
  return_op.reinit_domain_vector = reinit_vector;

  return_op.vmult = [](Range &v, const Range &u) {
    const auto mean = u.mean_value();

    v = u;
    v.add(-mean);
  };

  return_op.vmult_add = [](Range &v, const Range &u) {
    const auto mean = u.mean_value();

    v += u;
    v.add(-mean);
  };

  return_op.Tvmult     = return_op.vmult_add;
  return_op.Tvmult_add = return_op.vmult_add;

  return return_op;
}


/**
 * @relatesalso LinearOperator
 *
 * Return a LinearOperator that acts as a mean value filter. The vmult()
 * functions of this matrix subtract the mean values of the vector.
 *
 * The function takes a LinearOperator @p op and uses its range initializer
 * to create an mean value filter operator. The function also ensures that
 * the underlying Payload matches that of the input @p op.
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
mean_value_filter(const LinearOperator<Range, Domain, Payload> &op)
{
  auto return_op = mean_value_filter<Range, Payload>(op.reinit_range_vector);
  static_cast<Payload &>(return_op) = op.identity_payload();

  return return_op;
}


namespace internal
{
  namespace LinearOperatorImplementation
  {
    /**
     * A helper class that is responsible for the initialization of a vector
     * to be directly usable as the destination parameter, or source parameter
     * in an application of vmult of a matrix.
     *
     * The generic version of this class just calls
     * <code>Vector::reinit()</code> with the result of
     * <code>Matrix::m()</code> or <code>Matrix::n()</code>, respectively.
     * This class is specialized for more complicated data structures, such as
     * TrilinosWrappers::MPI::Vector, etc.
     */
    template <typename Vector>
    class ReinitHelper
    {
    public:
      /**
       * Initializes a vector v of the Range space to be directly usable as
       * the destination parameter in an application of vmult. Similar to the
       * reinit functions of the vector classes, the boolean determines
       * whether a fast initialization is done, i.e., if it is set to false the
       * content of the vector is set to 0.
       *
       * The generic version of this class just calls
       * <code>Vector::reinit()</code> with the result of
       * <code>Matrix::m()</code>.
       */
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &matrix,
                          Vector       &v,
                          bool          omit_zeroing_entries)
      {
        v.reinit(matrix.m(), omit_zeroing_entries);
      }

      /**
       * Initializes a vector of the Domain space to be directly usable as the
       * source parameter in an application of vmult. Similar to the reinit
       * functions of the vector classes, the boolean determines whether a
       * fast initialization is done, i.e., if it is set to false the content
       * of the vector is set to 0.
       *
       * The generic version of this class just calls
       * <code>Vector::reinit()</code> with the result of
       * <code>Matrix::n()</code>.
       */
      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &matrix,
                           Vector       &v,
                           bool          omit_zeroing_entries)
      {
        v.reinit(matrix.n(), omit_zeroing_entries);
      }
    };


    /**
     * A dummy class for LinearOperators that do not require any extensions to
     * facilitate the operations of the matrix.
     *
     * This is the Payload class typically associated with deal.II's native
     * SparseMatrix. To use Trilinos and PETSc sparse matrix classes it is
     * necessary to initialize a LinearOperator with their associated Payload.
     *
     *
     * @ingroup LAOperators
     */
    class EmptyPayload
    {
    public:
      /**
       * Default constructor
       *
       * Since this class does not do anything in particular and needs no
       * special configuration, we have only one generic constructor that can
       * be called under any conditions.
       */
      template <typename... Args>
      EmptyPayload(const Args &...)
      {}


      /**
       * Return a payload configured for identity operations
       */
      EmptyPayload
      identity_payload() const
      {
        return *this;
      }


      /**
       * Return a payload configured for null operations
       */
      EmptyPayload
      null_payload() const
      {
        return *this;
      }


      /**
       * Return a payload configured for transpose operations
       */
      EmptyPayload
      transpose_payload() const
      {
        return *this;
      }


      /**
       * Return a payload configured for inverse operations
       */
      template <typename Solver, typename Preconditioner>
      EmptyPayload
      inverse_payload(Solver &, const Preconditioner &) const
      {
        return *this;
      }
    };

    /**
     * Operator that returns a payload configured to support the addition of
     * two LinearOperators
     */
    inline EmptyPayload
    operator+(const EmptyPayload &, const EmptyPayload &)
    {
      return {};
    }

    /**
     * Operator that returns a payload configured to support the
     * multiplication of two LinearOperators
     */
    inline EmptyPayload
    operator*(const EmptyPayload &, const EmptyPayload &)
    {
      return {};
    }



    // A trait class that determines whether type T provides public
    // (templated or non-templated) vmult_add member functions
    template <typename Range, typename Domain, typename T>
    class has_vmult_add_and_Tvmult_add
    {
      template <typename C>
      static std::false_type
      test(...);

      template <typename C>
      static auto
      test(Range *r, Domain *d)
        -> decltype(std::declval<C>().vmult_add(*r, *d),
                    std::declval<C>().Tvmult_add(*d, *r),
                    std::true_type());

    public:
      // type is std::true_type if Matrix provides vmult_add and Tvmult_add,
      // otherwise it is std::false_type

      using type = decltype(test<T>(nullptr, nullptr));
    };


    // A helper function to apply a given vmult, or Tvmult to a vector with
    // intermediate storage
    template <typename Function, typename Range, typename Domain>
    void
    apply_with_intermediate_storage(Function      function,
                                    Range        &v,
                                    const Domain &u,
                                    bool          add)
    {
      GrowingVectorMemory<Range> vector_memory;

      typename VectorMemory<Range>::Pointer i(vector_memory);
      i->reinit(v, /*bool omit_zeroing_entries =*/true);

      function(*i, u);

      if (add)
        v += *i;
      else
        v = *i;
    }


    // A helper class to add a reduced matrix interface to a LinearOperator
    // (typically provided by Preconditioner classes)
    template <typename Range, typename Domain, typename Payload>
    class MatrixInterfaceWithoutVmultAdd
    {
    public:
      template <typename Matrix>
      void
      operator()(LinearOperator<Range, Domain, Payload> &op,
                 const Matrix                           &matrix)
      {
        op.vmult = [&matrix](Range &v, const Domain &u) {
          if (PointerComparison::equal(&v, &u))
            {
              // If v and u are the same memory location use intermediate
              // storage
              apply_with_intermediate_storage(
                [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
                v,
                u,
                /*bool add =*/false);
            }
          else
            {
              matrix.vmult(v, u);
            }
        };

        op.vmult_add = [&matrix](Range &v, const Domain &u) {
          // use intermediate storage to implement vmult_add with vmult
          apply_with_intermediate_storage(
            [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
            v,
            u,
            /*bool add =*/true);
        };

        op.Tvmult = [&matrix](Domain &v, const Range &u) {
          if (PointerComparison::equal(&v, &u))
            {
              // If v and u are the same memory location use intermediate
              // storage
              apply_with_intermediate_storage(
                [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
                v,
                u,
                /*bool add =*/false);
            }
          else
            {
              matrix.Tvmult(v, u);
            }
        };

        op.Tvmult_add = [&matrix](Domain &v, const Range &u) {
          // use intermediate storage to implement Tvmult_add with Tvmult
          apply_with_intermediate_storage(
            [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
            v,
            u,
            /*bool add =*/true);
        };
      }
    };


    // A helper class to add the full matrix interface to a LinearOperator
    template <typename Range, typename Domain, typename Payload>
    class MatrixInterfaceWithVmultAdd
    {
    public:
      template <typename Matrix>
      void
      operator()(LinearOperator<Range, Domain, Payload> &op,
                 const Matrix                           &matrix)
      {
        // As above ...

        MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>().operator()(
          op, matrix);

        // ... but add native vmult_add and Tvmult_add variants:

        op.vmult_add = [&matrix](Range &v, const Domain &u) {
          if (PointerComparison::equal(&v, &u))
            {
              apply_with_intermediate_storage(
                [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
                v,
                u,
                /*bool add =*/true);
            }
          else
            {
              matrix.vmult_add(v, u);
            }
        };

        op.Tvmult_add = [&matrix](Domain &v, const Range &u) {
          if (PointerComparison::equal(&v, &u))
            {
              apply_with_intermediate_storage(
                [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
                v,
                u,
                /*bool add =*/true);
            }
          else
            {
              matrix.Tvmult_add(v, u);
            }
        };
      }
    };
  } // namespace LinearOperatorImplementation
} // namespace internal


/**
 * @relatesalso LinearOperator
 *
 * A function that encapsulates generic @p matrix objects that act on a
 * compatible Vector type into a LinearOperator. The LinearOperator object
 * that is created stores a reference to the matrix object. Thus, @p matrix
 * must remain a valid reference for the whole lifetime of the LinearOperator
 * object.
 *
 * All changes made on @p matrix after the creation of the LinearOperator
 * object are reflected by the operator object. For example, it is a valid
 * procedure to first create a LinearOperator and resize, reassemble the
 * matrix later.
 *
 * The Matrix class in question must provide the following minimal interface:
 *
 * @code
 * class Matrix
 * {
 * public:
 *   // (type specific) information how to create a Range and Domain vector
 *   // with appropriate size and internal layout
 *
 *   // Application of matrix to vector src, writes the result into dst.
 *   vmult(Range &dst, const Domain &src);
 *
 *   // Application of the transpose of matrix to vector src, writes the
 *   // result into dst. (Depending on the usage of the linear operator
 *   // class this can be a dummy implementation throwing an error.)
 *   Tvmult(Range &dst, const Domain &src);
 * };
 * @endcode
 *
 * The following (optional) interface is used if available:
 *
 * @code
 * class Matrix
 * {
 * public:
 *   // Application of matrix to vector src, adds the result to dst.
 *   vmult_add(Range &dst, const Domain &src);
 *
 *   // Application of the transpose of matrix to vector src, adds the
 *   // result to dst.
 *   Tvmult_add(Range &dst, const Domain &src);
 * };
 * @endcode
 *
 * If the Matrix does not provide <code>vmult_add</code> and
 * <code>Tvmult_add</code>, they are implemented in terms of
 * <code>vmult</code> and <code>Tvmult</code> (requiring intermediate
 * storage).
 *
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload, typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const Matrix &matrix)
{
  // implement with the more generic variant below...
  return linear_operator<Range, Domain, Payload, Matrix, Matrix>(matrix,
                                                                 matrix);
}


/**
 * @relatesalso LinearOperator
 *
 * Variant of above function that takes an operator object @p
 * operator_exemplar as an additional reference. This object is used to
 * populate the reinit_domain_vector and reinit_range_vector function objects.
 * The reference @p matrix is used to construct vmult, Tvmult, etc.
 *
 * This variant can, for example, be used to encapsulate preconditioners (that
 * typically do not expose any information about the underlying matrix).
 *
 *
 * @ingroup LAOperators
 */
template <typename Range,
          typename Domain,
          typename Payload,
          typename OperatorExemplar,
          typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &operator_exemplar, const Matrix &matrix)
{
  using namespace internal::LinearOperatorImplementation;
  // Initialize the payload based on the input exemplar matrix
  LinearOperator<Range, Domain, Payload> return_op{
    Payload(operator_exemplar, matrix)};

  // Always store a reference to matrix and operator_exemplar in the lambda
  // functions. This ensures that a modification of the matrix after the
  // creation of a LinearOperator wrapper is respected - further a matrix
  // or an operator_exemplar cannot usually be copied...

  return_op.reinit_range_vector =
    [&operator_exemplar](Range &v, bool omit_zeroing_entries) {
      internal::LinearOperatorImplementation::ReinitHelper<
        Range>::reinit_range_vector(operator_exemplar, v, omit_zeroing_entries);
    };

  return_op.reinit_domain_vector = [&operator_exemplar](
                                     Domain &v, bool omit_zeroing_entries) {
    internal::LinearOperatorImplementation::ReinitHelper<
      Domain>::reinit_domain_vector(operator_exemplar, v, omit_zeroing_entries);
  };

  std::conditional_t<
    has_vmult_add_and_Tvmult_add<Range, Domain, Matrix>::type::value,
    MatrixInterfaceWithVmultAdd<Range, Domain, Payload>,
    MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>>()
    .
    operator()(return_op, matrix);

  return return_op;
}



/**
 * @relatesalso LinearOperator
 *
 * Variant of above function that takes a LinearOperator @p
 * operator_exemplar as an additional reference.
 * The reinit_domain_vector and reinit_range_vector function are copied
 * from the @p operator_exemplar object.
 *
 * The reference @p matrix is used to construct vmult, Tvmult, etc.
 *
 * This variant can, for example, be used to encapsulate preconditioners (that
 * typically do not expose any information about the underlying matrix).
 *
 *
 * @ingroup LAOperators
 */
template <typename Range, typename Domain, typename Payload, typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const LinearOperator<Range, Domain, Payload> &operator_exemplar,
                const Matrix                                 &matrix)
{
  using namespace internal::LinearOperatorImplementation;
  // Initialize the payload based on the LinearOperator exemplar
  auto return_op = operator_exemplar;

  std::conditional_t<
    has_vmult_add_and_Tvmult_add<Range, Domain, Matrix>::type::value,
    MatrixInterfaceWithVmultAdd<Range, Domain, Payload>,
    MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>>()
    .
    operator()(return_op, matrix);

  return return_op;
}


/** @} */

#ifndef DOXYGEN

//
// Ensure that we never capture a reference to a temporary by accident.
// to avoid "stack use after free".
//

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<Matrix>>>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &, Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<OperatorExemplar>>,
  typename = std::enable_if_t<
    !std::is_same_v<OperatorExemplar, LinearOperator<Range, Domain, Payload>>>>
LinearOperator<Range, Domain, Payload>
linear_operator(OperatorExemplar &&, const Matrix &) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<Matrix>>,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<OperatorExemplar>>,
  typename = std::enable_if_t<
    !std::is_same_v<OperatorExemplar, LinearOperator<Range, Domain, Payload>>>>
LinearOperator<Range, Domain, Payload>
linear_operator(OperatorExemplar &&, Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<Matrix>>>
LinearOperator<Range, Domain, Payload>
linear_operator(const LinearOperator<Range, Domain, Payload> &,
                Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<Matrix>>>
LinearOperator<Range, Domain, Payload>
linear_operator(Matrix &&) = delete;

template <
  typename Payload,
  typename Solver,
  typename Preconditioner,
  typename Range  = typename Solver::vector_type,
  typename Domain = Range,
  typename = std::enable_if_t<!std::is_lvalue_reference_v<Preconditioner>>,
  typename =
    std::enable_if_t<!std::is_same_v<Preconditioner, PreconditionIdentity>>,
  typename = std::enable_if_t<
    !std::is_same_v<Preconditioner, LinearOperator<Range, Domain, Payload>>>>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &,
                 Solver &,
                 Preconditioner &&) = delete;

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
