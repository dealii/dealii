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

#ifndef __deal2__linear_operator_h
#define __deal2__linear_operator_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector_memory.h>

#include <functional>
#include <type_traits>

#ifdef DEAL_II_WITH_CXX11

DEAL_II_NAMESPACE_OPEN

// Forward declarations:

template <typename Number> class Vector;

template <typename Range, typename Domain> class LinearOperator;
template <typename VECTOR = Vector<double>, typename MATRIX>
LinearOperator<VECTOR, VECTOR> linop(const MATRIX &matrix);

/**
 * A class to store the abstract concept of a linear operator.
 *
 * The primary purpose of this class is to provide syntactic sugar for
 * complex matrix-vector operations and free the user from having to
 * create, set up and handle intermediate storage locations by hand.
 *
 * This is achieved by storing the concept of such operations in an
 * "abstract" class LinearOperator that only holds knowledge on how to
 * apply a linear operation via an std::function object vmult, as well as,
 * information on how to create an appropriate Domain and Range vector.
 *
 * As an example consider the equation $(M+1/2\,k\,S)v=u$, where $M$
 * denotes the mass matrix and $S$ denotes the stiffness matrix. In order
 * to solve above equation with an iterative solver, one can write
 *
 * @code
 * dealii::SparseMatrix<double> mass_matrix;
 * dealii::SparseMatrix<double> stiffness_matrix;
 * const double k = ...;
 *
 * // Set up and assembly of mass and stiffness matrix
 *
 * const auto M = linop(mass_matrix);
 * const auto S = linop(stiffness_matrix);
 *
 * SolverControl solver_control (1000, 1e-12);
 * SolverCG<> solver (solver_control);
 *
 * Vector<double> v;
 * S.reinit_range_vector(v);
 * solver.solve(M + 0.5 * k * S, v, u, PreconditionIdentity());
 * @endcode
 *
 * This class is only available if deal.II was configured with C++11
 * support, i.e., if <code>DEAL_II_WITH_CXX11</code> is enabled during
 * cmake configure.
 *
 * @author: Luca Heltai, Matthias Maier, et al., 2015
 */
template <typename Range, typename Domain> class LinearOperator
{
public:

  /**
   * Create an empty LinearOperator object
   */
  LinearOperator ()
  {
    vmult = [](Range &, const Domain &)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::vmult called"));
    };

    vmult_add = [](Range &, const Domain &)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::vmult_add called"));
    };

    Tvmult = [](Domain &, const Range &)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::Tvmult called"));
    };

    Tvmult_add = [](Domain &, const Range &)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::Tvmult_add called"));
    };

    reinit_range_vector = [](Range &, bool)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::reinit_range_vector method called"));
    };

    reinit_domain_vector = [](Domain &, bool)
    {
      Assert(false, ExcMessage("Uninitialized LinearOperator<Range, "
                               "Domain>::reinit_domain_vector method called"));
    };
  }

  /**
   * Create a LinearOperator object from an object @p op for which a
   * conversion function <code>linop</code> is defined.
   */
  template<typename Op>
  LinearOperator (const Op &op)
  {
    *this = linop<Range>(op);
  }

  /**
   * Copy assignment operator for an object @p op for which a conversion
   * function <code>linop</code> is defined.
   */
  template <typename Op>
  LinearOperator<Range, Domain> &operator=(const Op &op)
  {
    return *this = linop<Range>(op);
  }

  /**
   * Application of the LinearOperator object to a vector @p u of the
   * Domain space giving a vector @p v of the Range space.
   */
  std::function<const void(Range &v, const Domain &u)> vmult;

  /**
   * Application of the LinearOperator object to a vector @p u of the
   * Domain space. The result is added to the vector @p v.
   */
  std::function<const void(Range &v, const Domain &u)> vmult_add;

  /**
   * Application of the transpose LinearOperator object to a vector @p u of
   * the Range space giving a vector @p v of the Domain space.
   */
  std::function<const void(Domain &v, const Range &u)> Tvmult;

  /**
   * Application of the transpose LinearOperator object to a vector @p u of
   * the Range space. The result is added to the vector @p v.
   */
  std::function<const void(Domain &v, const Range &u)> Tvmult_add;

  /**
   * Initializes a vector of the Range space to be directly usable as the
   * destination parameter in an application of vmult. This function object
   * does not set the content of the vector to 0.
   */
  std::function<void(Range &, bool)> reinit_range_vector;

  /**
   * Initializes a vector of the Range space to be directly usable as the
   * source parameter in an application of vmult. This function object does
   * not set the content of the vector to 0.
   */
  std::function<void(Domain &, bool)> reinit_domain_vector;

  /**
   * A memory object for vectors of Range space used for intermediate
   * storage during computations in vmult.
   */
  mutable GrowingVectorMemory<Range> range_vector_memory;

  /**
   * A memory object for vectors of Domain space used for intermediate
   * storage during computations in vmult.
   */
  mutable GrowingVectorMemory<Domain> domain_vector_memory;


  /**
   * Addition with a LinearOperator @p second_op with the same Domain
   * and Range.
   */
  LinearOperator<Range, Domain> &
  operator+=(const LinearOperator<Range, Domain> &second_op)
  {
    *this = *this + second_op;
    return *this;
  }

  /**
   * Subtraction with a LinearOperator @p second_op with the same Domain
   * and Range.
   */
  LinearOperator<Range, Domain> &
  operator-=(const LinearOperator<Range, Domain> &second_op)
  {
    *this = *this - second_op;
    return *this;
  }

  /**
   * Concatenation of the LinearOperator with an endormorphism @p second_op
   * on the Domain space.
   */
  LinearOperator<Range, Domain> &
  operator*=(const LinearOperator<Domain, Domain> &second_op)
  {
    *this = *this * second_op;
    return *this;
  }
};



/**
 * Addition of two linear operators @p first_op and @p second_op given
 * by $(\text{first_op}+\text{second_op})x:=\text{first_op}(x)+\text{second_op}(x)$
 */
template <typename Range, typename Domain>
LinearOperator<Range, Domain>
operator+(const LinearOperator<Range, Domain> &first_op,
          const LinearOperator<Range, Domain> &second_op)
{
  LinearOperator<Range, Domain> return_op;

  return_op.reinit_range_vector = first_op.reinit_range_vector;
  return_op.reinit_domain_vector = first_op.reinit_domain_vector;

  // ensure to have valid computation objects by catching first_op and
  // second_op by value

  return_op.vmult = [first_op, second_op](Range &v, const Domain &u)
  {
    first_op.vmult(v, u);
    second_op.vmult_add(v, u);
  };

  return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u)
  {
    first_op.vmult_add(v, u);
    second_op.vmult_add(v, u);
  };

  return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u)
  {
    second_op.Tvmult(v, u);
    first_op.Tvmult_add(v, u);
  };

  return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u)
  {
    second_op.Tvmult_add(v, u);
    first_op.Tvmult_add(v, u);
  };

  return return_op;
}


/**
 * Subtraction of two linear operators @p first_op and @p second_op given
 * by $(\text{first_op}-\text{second_op})x:=\text{first_op}(x)-\text{second_op}(x)$
 */
template <typename Range, typename Domain>
LinearOperator<Range, Domain>
operator-(const LinearOperator<Range, Domain> &first_op,
          const LinearOperator<Range, Domain> &second_op)
{
  // implement with addition and scalar multiplication
  return first_op + (-1. * second_op);
}


/**
 * Concatenation of two linear operators @p first_op and @p second_op given
 * by $(\text{first_op}*\text{second_op})x:=\text{first_op}(\text{second_op}(x))$
 */
template <typename Range, typename Intermediate, typename Domain>
LinearOperator<Range, Domain>
operator*(const LinearOperator<Range, Intermediate> &first_op,
          const LinearOperator<Intermediate, Domain> &second_op)
{
  LinearOperator<Range, Domain> return_op;

  return_op.reinit_domain_vector = second_op.reinit_domain_vector;
  return_op.reinit_range_vector = first_op.reinit_range_vector;

  // ensure to have valid computation objects by catching first_op and
  // second_op by value

  // Note: The range_vector_memory and domain_vector_memory objects inside
  // the lambda exrepessions are the one belonging to the by-value captured
  // objects first_op and second_op. Thus, they belong in essence to the
  // std::function object vmult (vmult_add, ...) - and therefore have the
  // same lifetime as the function object and not the LinearOperator
  // parameters first_op and second_op.

  return_op.vmult = [first_op, second_op](Range &v, const Domain &u)
  {
    Intermediate *i = second_op.range_vector_memory.alloc();
    second_op.reinit_range_vector(*i, /*bool fast =*/ true);
    second_op.vmult(*i, u);
    first_op.vmult(v, *i);
    second_op.range_vector_memory.free(i);
  };

  return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u)
  {
    Intermediate *i = second_op.range_vector_memory.alloc();
    second_op.reinit_range_vector(*i, /*bool fast =*/ true);
    second_op.vmult(*i, u);
    first_op.vmult_add(v, *i);
    second_op.range_vector_memory.free(i);
  };

  return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u)
  {
    Intermediate *i = first_op.domain_vector_memory.alloc();
    first_op.reinit_domain_vector(*i, /*bool fast =*/ true);
    first_op.Tvmult(*i, u);
    second_op.Tvmult(v, *i);
    first_op.domain_vector_memory.free(i);
  };

  return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u)
  {
    Intermediate *i = first_op.domain_vector_memory.alloc();
    first_op.reinit_domain_vector(*i, /*bool fast =*/ true);
    first_op.Tvmult(*i, u);
    second_op.Tvmult_add(v, *i);
    first_op.domain_vector_memory.free(i);
  };

  return return_op;
}


/**
 * Scalar multiplication of a ScalarOperator object from the left.
 *
 * The Domain and Range types must implement the following
 * <code>operator*=</code> member functions accepting the appropriate
 * scalar Number type for rescaling:
 *
 * @code
 * Domain & operator *=(Number);
 * Range & operator *=(Number);
 * @endcode
 */
template <typename Range, typename Domain>
LinearOperator<Range, Domain>
operator*(typename Range::value_type  number,
          const LinearOperator<Range, Domain> &op)
{
  static_assert(
    std::is_convertible<typename Range::value_type, typename Domain::value_type>::value,
    "Range and Domain must have implicitly convertible 'value_type's");

  LinearOperator<Range, Domain> return_op = op;

  // ensure to have valid computation objects by catching number and op by
  // value

  return_op.vmult = [number, op](Range &v, const Domain &u)
  {
    op.vmult(v,u);
    v *= number;
  };

  return_op.vmult_add = [number, op](Range &v, const Domain &u)
  {
    v /= number;
    op.vmult_add(v,u);
    v *= number;
  };

  return_op.Tvmult = [number, op](Domain &v, const Range &u)
  {
    op.Tvmult(v,u);
    v *= number;
  };

  return_op.Tvmult_add = [number, op](Domain &v, const Range &u)
  {
    v /= number;
    op.Tvmult_add(v,u);
    v *= number;
  };

  return return_op;
}


/**
 * Scalar multiplication of a ScalarOperator object from the right.
 *
 * The Domain and Range types must implement the following
 * <code>operator*=</code> member functions for rescaling:
 *
 * @code
 * Domain & operator *=(Number);
 * Range & operator *=(Number);
 * @endcode
 */
template <typename Range, typename Domain>
LinearOperator<Range, Domain>
operator*(const LinearOperator<Range, Domain> &op,
          typename Domain::value_type  number)
{
  static_assert(
    std::is_convertible<typename Range::value_type, typename Domain::value_type>::value,
    "Range and Domain must have implicitly convertible 'value_type's");

  return number * op;
}


/**
 * Returns the transpose linear operations of @p op.
 */
template <typename Range, typename Domain>
LinearOperator<Domain, Range>
transpose_linop(const LinearOperator<Range, Domain> &op)
{
  LinearOperator<Domain, Range> return_op;

  return_op.reinit_range_vector = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = op.Tvmult;
  return_op.vmult_add = op.Tvmult_add;
  return_op.Tvmult = op.vmult;
  return_op.Tvmult_add = op.vmult_add;

  return return_op;
}


/**
 * Returns a LinearOperator that is the identity of the vector space
 * Domain.
 *
 * The function takes an <code>std::function</code> object @p
 * exemplar as an argument to initialize the reinit_range_vector
 * and reinit_domain_vector objects of the LinearOperator object.
 */
template <typename Range>
LinearOperator<Range, Range>
identity_linop(const std::function<void(Range &, bool)> &exemplar)
{
  LinearOperator<Range, Range> return_op;

  return_op.reinit_range_vector = exemplar;
  return_op.reinit_domain_vector = exemplar;

  return_op.vmult = [](Range &v, const Range &u)
  {
    Assert(u.size() == v.size(), ExcDimensionMismatch(u.size(), v.size()));
    v = u;
  };

  return_op.vmult_add = [](Range &v, const Range &u)
  {
    Assert(u.size() == v.size(), ExcDimensionMismatch(u.size(), v.size()));
    v += u;
  };

  return_op.Tvmult = [](Range &v, const Range &u)
  {
    Assert(u.size() == v.size(), ExcDimensionMismatch(u.size(), v.size()));
    v = u;
  };

  return_op.Tvmult_add = [](Range &v, const Range &u)
  {
    Assert(u.size() == v.size(), ExcDimensionMismatch(u.size(), v.size()));
    v += u;
  };

  return return_op;
}


/**
 * Returns an object representing the inverse of the LinearOperator @p op.
 *
 * The function takes references @p solver and @p preconditioner to an
 * iterative solver and a preconditioner that are used in the
 * <code>vmult</code> and <code>Tvmult</code> implementations of the
 * LinearOperator object.
 *
 * The LinearOperator object that is created stores a reference to @p
 * solver and @p preconditioner. Thus, both objects must remain a valid
 * reference for the whole lifetime of the LinearOperator object. Internal
 * data structures of the @p solver object will be modified upon invocation
 * of <code>vmult</code> or <code>Tvmult</code>.
 */
template <typename Solver, typename Preconditioner>
LinearOperator<typename Solver::vector_type, typename Solver::vector_type>
inverse_linop(Solver &solver,
              const Preconditioner &preconditioner,
              const LinearOperator<typename Solver::vector_type, typename Solver::vector_type> &op)
{
  typedef typename Solver::vector_type Vector;

  LinearOperator<Vector, Vector> return_op;

  return_op.reinit_range_vector = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = [op, &solver, &preconditioner](Vector &v, const Vector &u)
  {
    solver.solve(op, v, u, preconditioner);
  };

  // Note: The range_vector_memory and domain_vector_memory objects inside
  // the lambda exrepessions are the one belonging to the by-value captured
  // objects first_op and second_op. Thus, they belong in essence to the
  // std::function object vmult (vmult_add, ...) - and therefore have the
  // same lifetime as the function object and not the LinearOperator
  // parameter op.

  return_op.vmult_add =
    [op, &solver, &preconditioner](Vector &v, const Vector &u)
  {
    Vector *v2 = op.range_vector_memory.alloc();
    op.reinit_range_vector(*v2, /*bool fast =*/ true);
    solver.solve(op, *v2, u, preconditioner);
    v += *v2;
    op.range_vector_memory.free(v2);
  };

  return_op.Tvmult = [op, &solver, &preconditioner]( Vector &v, const
                                                     Vector &u)
  {
    solver.solve(transpose_linop(op), v, u, preconditioner);
  };

  return_op.Tvmult =
    [op, &solver, &preconditioner](Vector &v, const Vector &u)
  {
    Vector *v2 = op.range_vector_memory.alloc();
    op.reinit_range_vector(*v2, /*bool fast =*/ true);
    solver.solve(transpose_linop(op), *v2, u, preconditioner);
    v += *v2;
    op.range_vector_memory.free(v2);
  };

  return return_op;
}


/**
 * A factory class that is responsible to create a reinit_range_vector
 * object for a given pair of vector type Range and matrix type Matrix.
 *
 * The generic version of this class just calls Range::reinit() with the
 * result of Matrix::m(). This class is specialized for more complicated
 * data structures, such as TrilinosWrappers::MPI::Vector, etc.
 */
template<typename Range>
class ReinitRangeFactory
{
public:

  template <typename Matrix>
  std::function<void(Range &, bool)>
  operator()(const Matrix &matrix)
  {
    return [&matrix](Range &v, bool fast)
    {
      v.reinit(matrix.m(), fast);
    };
  }
}


/**
 * A factory class that is responsible to create a reinit_domain_vector
 * object for a given pair of vector type Domain and matrix type Matrix.
 *
 * The generic version of this class just calls Domain::reinit() with the
 * result of Matrix::m(). This class is specialized for more complicated
 * data structures, such as TrilinosWrappers::MPI::Vector, etc.
 */
template<typename Domain>
class ReinitDomainFactory
{
public:

  template <typename Matrix>
  std::function<void(Domain &, bool)>
  operator()(const Matrix &matrix)
  {
    return [&matrix](Domain &v, bool fast)
    {
      v.reinit(matrix.m(), fast);
    };
  }
}


/**
 * A function that encapsulates generic @p matrix objects that act on a
 * compatible Vector type  into a LinearOperator. The LinearOperator object
 * that is created stores a reference to the matrix object. Thus, @p matrix
 * must remain a valid reference for the whole lifetime of the
 * LinearOperator object.
 *
 * All changes made on @p matrix after the creation of the LinearOperator
 * object are reflected by the operator object. For example, it is a valid
 * procedure to first create a LinearOperator and resize, reassemble the
 * matrix later.
 *
 * The Matrix class in question must provide the following minimal
 * interface:
 *
 * @code
 * class Matrix
 * {
 * public:
 *
 *   // (Type specific) information how to create a VECTOR with appropriate
 *   // size (or internal layout) for domain and range
 *
 *   // Application of matrix to vector src. Writes result into dst
 *   vmult(VECTOR &dst, const VECTOR &src);
 *
 *   // Application of matrix to vector src. Add result to dst
 *   vmult_add(VECTOR &dst, const VECTOR &src);
 *
 *   // Application of transpose of matrix to vector src. Writes result
 *   // into dst. (Depending on the usage of the linear operator class this
 *   // can be a dummy implementation throwing an error)
 *   Tvmult(VECTOR &dst, const VECTOR &src);
 *
 *   // Application of matrix to vector src. Add result to dst
 *   // Application of transpose of matrix to vector src. Add result to
 *   /  dst. (Depending on the usage of the linear operator class this
 *   // can be a dummy implementation throwing an error)
 *   Tvmult_add(VECTOR &dst, const VECTOR &src);
 * };
 * @endcode
 *
 * @author: Matthias Maier, 2015
 */
template <typename Range = Vector<double>, typename Domain = Range,
          typename Matrix>
LinearOperator<Range, Domain> linop(const Matrix &Matrix)
{
  LinearOperator<Range, Domain> return_op;

  // Store all information about Matrix in the lambda functions by
  // reference to ensure that a modification of Matrix after the creation
  // of a LinearOperator wrapper is respected.

  return_op.reinit_range_vector =
      ReinitRangeFactory<Range>.operator()(matrix);

  return_op.reinit_domain_vector =
      ReinitDomainFactory<Domain>.operator()(matrix);

  return_op.vmult = [&matrix](Range &v, const Domain &u)
  {
    matrix.vmult(v,u);
  };

  return_op.vmult_add = [&matrix](Range &v, const Domain &u)
  {
    matrix.vmult_add(v,u);
  };

  return_op.Tvmult = [&matrix](Domain &v, const Range &u)
  {
    matrix.Tvmult(v,u);
  };

  return_op.Tvmult_add = [&matrix](Domain &v, const Range &u)
  {
    matrix.Tvmult_add(v,u);
  };

  return return_op;
}


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_CXX11
#endif
