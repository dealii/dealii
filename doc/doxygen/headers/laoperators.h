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
 * <code>DEAL_II_WITH_CXX11=on</code> during configuration) a versatile
 * mechanism for storing the concept of a linear operator is available.
 *
 * This is done with a LinearOperator class that, similarly to the abstract
 * MATRIX interface, defines a minimal interface for <i>applying</i> a
 * linear operation on a vector.
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
 * @ingroup LAC
 * @ingroup MATRICES
 */
