// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



/**
 * @defgroup Preconditioners Preconditioners and Relaxation Operators
 *
 * <h3>Preconditioners</h3>
 *
 * Preconditioners are used to accelerate the iterative solution of linear
 * systems. Typical preconditioners are Jacobi, Gauss-Seidel, or SSOR, but the
 * library also supports more complex ones such as Vanka or incomplete LU
 * decompositions (ILU). In addition, sparse direct solvers can be used as
 * preconditioners when available.
 *
 * Broadly speaking, preconditioners are operators, which are multiplied with
 * a matrix to improve conditioning. The idea is, that the preconditioned
 * system <i>P<sup>-1</sup>Ax = P<sup>-1</sup>b</i> is much easier to solve
 * than the original system <i>Ax = b</i>. What this means exactly depends on
 * the structure of the matrix and cannot be discussed here in generality. For
 * symmetric, positive definite matrices <i>A</i> and <i>P</i>, it means that
 * the spectral condition number (the quotient of greatest and smallest
 * eigenvalue) of <i>P<sup>-1</sup>A</i> is much smaller than the one of
 * <i>A</i>.
 *
 * At hand of the simplest example, Richardson iteration, implemented in
 * SolverRichardson, the preconditioned iteration looks like
 * @f[
 *  x^{k+1} = x^k - P^{-1} \bigl(A x^k - b\bigr).
 * @f]
 * Accordingly, preconditioning amounts to applying a linear operator to the
 * residual, and consequently, the action of the preconditioner
 * <i>P<sup>-1</sup></i> is implemented as <tt>vmult()</tt>.
 * Templates in deal.II that require a preconditioner indicate the
 * requirement with
 * @ref ConceptPreconditionerType "the PreconditionerType concept". In
 * practice, one can usually treat any matrix-like object which defines
 * <code>vmult()</code> and <code>Tvmult()</code> as a preconditioner. All
 * preconditioner classes in this group implement this interface.
 *
 * When used
 * in Krylov space methods, it is up to the method, whether it simply
 * replaces multiplications with <i>A</i> by those with
 * <i>P<sup>-1</sup>A</i> (for instance SolverBicgstab), or does more
 * sophisticated things. SolverCG for instance uses
 * <i>P<sup>-1</sup></i> to define an inner product, which is the
 * reason why it requires a symmetric, positive definite operator <i>P</i>.
 *
 * <h3>Relaxation methods</h3>
 *
 * Many preconditioners rely on an additive splitting <i>A = P - N</i>
 * into two matrices. In this case, the iteration step of the
 * Richardson method above can be simplified to
 * @f[
 *  x^{k+1} = P^{-1} \bigl(N x^k + b\bigr),
 * @f]
 * thus avoiding multiplication with <i>A</i> completely. We call
 * operators mapping the previous iterate <i>x<sup>k</sup></i> to the
 * next iterate in this way relaxation operators. Their generic
 * interface is given by @ref ConceptRelaxationType "the RelaxationType concept".
 * The classes with names starting with <tt>Relaxation</tt> in this group
 * implement this interface, as well as the preconditioners
 * PreconditionJacobi, PreconditionSOR, PreconditionBlockJacobi,
 * PreconditionBlockSOR, and PreconditionBlockSSOR.
 *
 * <h3>The interface</h3>
 *
 * In this section, we discuss the interface preconditioners usually
 * have to provide to work inside the deal.II library.
 *
 * <h4>Initialization</h4>
 *
 * In order to be able to be stored in containers, all preconditioners
 * have a constructor with no arguments. Since this will typically
 * produce a useless object, all preconditioners have a function
 * @code
 *   void initialize (...)
 * @endcode
 *
 * This function receives the matrix to be preconditioned as well as
 * additional required parameters and sets up the internal structures
 * of the preconditioner.
 *
 * <h4>Relaxation methods</h4>
 *
 * Some preconditioners, like SOR and Jacobi, were used as iterative solvers
 * long before they were used as preconditioners. Thus, they satisfy both
 * @ref ConceptMatrixType "MatrixType" and
 * @ref ConceptRelaxationType "RelaxationType" concepts.
 *
 * @ingroup LAC
 * @ingroup Matrices
 */
