// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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



/**
 * @defgroup LAOperators Linear Operators
 *
 * <h3>Linear Operator</h3>
 *
 * If deal.II is configured with C++11 support (i.e.,
 * <code>DEAL_II_WITH_CXX11=ON</code> or <code>DEAL_II_WITH_CXX14=ON</code>
 * during configuration) a versatile mechanism for storing the concept of a
 * linear operator is available. (For questions about C++11, see
 * @ref CPP11 .)
 *
 * This is done with a LinearOperator class that, like
 * @ref ConceptMatrixType "the MatrixType concept",
 * defines a minimal interface for <i>applying</i> a linear operation on a
 * vector.
 *
 * @code
 *   std::function<void(Range &, const Domain &)> vmult;
 *   std::function<void(Range &, const Domain &)> vmult_add;
 *   std::function<void(Domain &, const Range &)> Tvmult;
 *   std::function<void(Domain &, const Range &)> Tvmult_add;
 * @endcode
 *
 * Thus, such an object can be used as a matrix object in all
 * @ref Solvers "iterative solver" classes, either as a matrix object, or as
 * @ref Preconditioners "preconditioner".
 *
 * The big advantage of the LinearOperator class is that it provides
 * syntactic sugar for complex matrix-vector operations. As an example
 * consider the operation $(A+k\,B)\,C$, where $A$, $B$ and $C$ denote
 * (possibly different) SparseMatrix objects. In order to construct a
 * LinearOperator <code>op</code> that performs above computation when
 * applied on a vector, one can write:
 * @code
 * dealii::SparseMatrix<double> A, B, C;
 * double k;
 * // Setup and assembly...
 *
 * const auto op_a = linear_operator(A);
 * const auto op_b = linear_operator(B);
 * const auto op_c = linear_operator(C);
 *
 * const auto op = (op_a + k * op_b) * op_c;
 * @endcode
 * Now, <code>op</code> can be used as a matrix object for further
 * computation.
 *
 * The linear_operator() function can be used to wrap an ordinary matrix or
 * preconditioner object into a LinearOperator. A linear operator can be
 * transposed with transpose_operator(), or inverted by using the
 * inverse_operator() together with an iterative solver.
 *
 * For objects of type LinearOperator, all vector space operations, i.e.,
 * addition and subtraction, scalar multiplication and composition (of
 * compatible linear operators) are implemented:
 * @code
 * dealii::LinearOperator<> op_a, op_b;
 * double k;
 *
 * // vector space addition, subtraction and scalar multiplication
 * op_a + op_b;
 * op_a - op_b;
 * k * op_a;
 * op_a * k;
 *
 * // in-place variants
 * op_a += op_b;
 * op_a -= op_b;
 * op_a *= k;
 *
 * // operator composition
 * op_a * op_b;
 * op_a *= op_b; // If op_b is an endomorphism of the domain space of op_a
 * @endcode
 *
 * block_operator() and block_diagonal_operator() provide further
 * encapsulation of individual linear operators into blocked linear
 * operator variants.
 *
 * @note The LinearOperator facility obsoletes some of the @ref Matrix2
 * "derived matrix" classes, such as BlockDiagonalMatrix, IterativeInverse,
 * ProductMatrix, ScaledMatrix, ProductSparseMatrix,
 * InverseMatrixRichardson, SchurMatrix, ShiftedMatrix,
 * ShiftedMatrixGeneralized, TransposeMatrix
 *
 * @note As explained below, when using LinearOperator as <code>res = op_a*x</code>
 * a PackagedOperation class instance is generated behind-the-curtains.
 * Consequently, the user program has to include header files for both classes
 * for compilation to be successful. In an attempt to make easier the decision of which
 * headers to include in what circumstances and to prevent hidden templates-related
 * compiler errors, all headers relevant to LinearOperator are grouped in
 * <deal.ii/lac/linear_operator_tools.h>.
 *
 * <h3>Packaged Operation</h3>
 *
 * An  application of a LinearOperator object to a vector via
 * <code>operator*</code> yields a PackagedOperation object that stores
 * this computation.
 *
 * The PackagedOperation class allows lazy evaluation of expressions
 * involving vectors and linear operators. This is done by storing the
 * computational expression and only performing the computation when either
 * the object is implicitly converted to a vector object, or
 * PackagedOperation::apply() (or PackagedOperation::apply_add()) is
 * invoked by hand. This avoids unnecessary temporary storage of
 * intermediate results.
 *
 * As an example consider the addition of multiple vectors:
 * @code
 *   dealii::Vector<double> a, b, c, d;
 *   // ..
 *   dealii::Vector<double> result = a + b - c + d;
 * @endcode
 * Converting the PackagedOperation <code>a + b - c + d</code> to a vector
 * results in code equivalent to the following code
 * @code
 *   dealii::Vector<double> a, b, c, d;
 *   // ..
 *   dealii::Vector<double> result = a;
 *   result += b;
 *   result -= c;
 *   result += d;
 * @endcode
 * that avoids any intermediate storage. As a second example (involving a
 * LinearOperator object) consider the computation of a residual $b-Ax$:
 *
 * @code
 *   dealii::SparseMatrix<double> A;
 *   dealii::Vector<double> b, x;
 *   // ..
 *   const auto op_a = linear_operator(A);
 *
 *   dealii::Vector<double> residual =  b - op_a * x;
 * @endcode
 * Here, the expression <code>b - op_a * x</code> results again in an
 * object of type PackagedOperation that stores the <i>sequence of
 * operations</i> that should be performed using the two vectors and the
 * linear operator. Converting the expression to a vector (as happens here
 * with the assignment to the vector <code>residual</code>) executes the
 * computation (see the following note).
 *
 * @note
 * Lazy evaluation of a computational expression necessarily involves
 * references to the underlying vector and matrix objects. For example, the
 * creation of a <code>residual_expr</code> object
 * @code
 *   auto residual_expr =  b - op_a * x;
 * @endcode
 * stores the computational expression of the residual with references to
 * the vector <code>b</code> and matrix <code>A</code>. It does not perform
 * any computation at this point. In particular, if <code>b</code> or
 * <code>A</code> are changed <b>after</b> the creation of
 * <code>residual_expr</code> every subsequent evaluation of the expression
 * is performed with the new values
 * @code
 *   auto residual_expr =  b - op_a * x;
 *   residual_expr.apply(tmp);  // tmp is a Vector<double>
 *
 *   // modify b, or A
 *
 *   residual_expr.apply(tmp2); // tmp2 is a Vector<double>
 *
 *   // tmp and tmp2 are different
 * @endcode
 * Thus, as a safeguard, if you want to compute the result of an expression
 * right away, always explicitly use a vector type on the left side (and
 * not <code>auto</code>):
 * @code
 *   Vector<double> residual =  b - op_a * x; // computes the residual at this point
 * @endcode
 *
 *
 * @ingroup LAC
 * @ingroup MATRICES
 */
