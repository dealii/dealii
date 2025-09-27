// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_schur_complement_h
#define dealii_schur_complement_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/vector_memory.h>


DEAL_II_NAMESPACE_OPEN

/**
 * @name Creation of a LinearOperator related to the Schur Complement
 */
/** @{ */

/**
 * @relatesalso LinearOperator
 *
 * Return a LinearOperator that performs the operations associated with the
 * Schur complement. There are two additional helper functions,
 * condense_schur_rhs() and postprocess_schur_solution(), that are likely
 * necessary to be used in order to perform any useful tasks in linear algebra
 * with this operator.
 *
 * We construct the definition of the Schur complement in the following way:
 *
 * Consider a general system of linear equations that can be decomposed into
 * two major sets of equations:
 * @f{eqnarray*}{
 * \mathbf{K}\mathbf{d} = \mathbf{f}
 * \quad \Rightarrow\quad
 * \left(\begin{array}{cc}
 *    A & B \\ C & D
 * \end{array}\right)
 * \left(\begin{array}{cc}
 *    x \\ y
 * \end{array}\right)
 * =
 * \left(\begin{array}{cc}
 *    f \\ g
 * \end{array}\right),
 * @f}
 * where $ A,B,C,D $  represent general subblocks of the matrix $ \mathbf{K} $
 * and, similarly, general subvectors of $ \mathbf{d},\mathbf{f} $ are given
 * by $ x,y,f,g $ .
 *
 * This is equivalent to the following two statements:
 * @f{eqnarray*}{
 *   (1) \quad Ax + By &=& f \\
 *   (2) \quad Cx + Dy &=& g \quad .
 * @f}
 *
 * Assuming that $ A,D $ are both square and invertible, we could then perform
 * one of two possible substitutions,
 * @f{eqnarray*}{
 *   (3) \quad x &=& A^{-1}(f - By) \quad \text{from} \quad (1) \\
 *   (4) \quad y &=& D^{-1}(g - Cx) \quad \text{from} \quad (2) ,
 * @f}
 * which amount to performing block Gaussian elimination on this system of
 * equations.
 *
 * For the purpose of the current implementation, we choose to substitute (3)
 * into (2)
 * @f{eqnarray*}{
 *   C \: A^{-1}(f - By) + Dy &=& g \\
 *   -C \: A^{-1} \: By + Dy &=& g - C \: A^{-1} \: f \quad .
 * @f}
 * This leads to the result
 * @f[
 *   (5) \quad (D - C\: A^{-1} \:B)y  = g - C \: A^{-1} f
 *       \quad \Rightarrow \quad Sy = g'
 * @f]
 * with $ S = (D - C\: A^{-1} \:B) $ being the Schur complement and the
 * modified right-hand side vector $ g' = g - C \: A^{-1} f $ arising from the
 * condensation step. Note that for this choice of $ S $, submatrix $ D $ need
 * not be invertible and may thus be the null matrix. Ideally $ A $ should be
 * well-conditioned.
 *
 * So for any arbitrary vector $ a $, the Schur complement performs the
 * following operation:
 * @f[
 *   (6) \quad Sa = (D - C \: A^{-1} \: B)a
 * @f]
 *
 * A typical set of steps needed the solve a linear system (1),(2) would be:
 * 1. Define the inverse matrix @p A_inv (using inverse_operator()).
 * 2. Define the Schur complement $ S $ (using schur_complement()).
 * 3. Define iterative inverse matrix $ S^{-1} $ such that (6)
 * holds. It is necessary to use a solver with a preconditioner to compute the
 * approximate inverse operation of $ S $ since we never compute $ S $
 * directly, but rather the result of its operation. To achieve this, one may
 * again use the inverse_operator() in conjunction with the Schur complement
 * that we've just constructed. Observe that the both $ S $ and its
 * preconditioner operate over the same space as $ D $.
 * 4. Perform pre-processing step on the RHS of (5) using
 * condense_schur_rhs():
 *    @f[
 *      g' = g - C \: A^{-1} \: f
 *    @f]
 * 5. Solve for $ y $ in (5):
 *    @f[
 *      y =  S^{-1} g'
 *    @f]
 * 6. Perform the post-processing step from (3) using
 * postprocess_schur_solution():
 *    @f[
 *      x =  A^{-1} (f - By)
 *    @f]
 *
 * An illustration of typical usage of this operator for a fully coupled
 * system is given below.
 * @code
 * #include <deal.II/lac/schur_complement.h>
 *
 * // Given BlockMatrix K and BlockVectors d,F
 *
 * // Decomposition of tangent matrix
 * const auto A = linear_operator(K.block(0,0));
 * const auto B = linear_operator(K.block(0,1));
 * const auto C = linear_operator(K.block(1,0));
 * const auto D = linear_operator(K.block(1,1));
 *
 * // Decomposition of solution vector
 * auto x = d.block(0);
 * auto y = d.block(1);
 *
 * // Decomposition of RHS vector
 * auto f = F.block(0);
 * auto g = F.block(1);
 *
 * // Construction of inverse of Schur complement
 * const auto prec_A = PreconditionSelector<...>(A);
 * const auto A_inv = inverse_operator<...>(A,prec_A);
 * const auto S = schur_complement(A_inv,B,C,D);
 *
 * // D and S operate on same space
 * const auto S_prec = PreconditionSelector<...>(D);
 * const auto S_inv = inverse_operator<...>(S,...,prec_S);
 *
 * // Solve reduced block system
 * // PackagedOperation that represents the condensed form of g
 * auto rhs = condense_schur_rhs (A_inv,C,f,g);
 *
 * // Solve for y
 * y = S_inv * rhs;
 *
 * // Compute x using resolved solution y
 * x = postprocess_schur_solution (A_inv,B,y,f);
 * @endcode
 *
 * In the above example, the preconditioner for $ S $ was defined as the
 * preconditioner for $ D $, which is valid since they operate on the same
 * space. However, if $ D $ and $ S $ are too dissimilar, then this may lead
 * to a large number of solver iterations as $ \text{prec}(D) $ is not a good
 * approximation for $ S^{-1} $.
 *
 * A better preconditioner in such a case would be one that provides a more
 * representative approximation for $ S^{-1} $. One approach is shown in
 * step-22, where $ D $ is the null matrix and the preconditioner for $ S^{-1}
 * $ is derived from the @ref GlossMassMatrix "mass matrix" over this space.
 *
 * From another viewpoint, a similar result can be achieved by first
 * constructing an object that represents an approximation for $ S $ wherein
 * expensive operation, namely $ A^{-1} $, is approximated. Thereafter we
 * construct the approximate inverse operator $ \tilde{S}^{-1} $ which is then
 * used as the preconditioner for computing $ S^{-1} $.
 * @code
 * // Construction of approximate inverse of Schur complement
 * const auto A_inv_approx = linear_operator(preconditioner_A);
 * const auto S_approx = schur_complement(A_inv_approx,B,C,D);
 *
 * // D and S_approx operate on same space
 * const auto S_approx_prec = PreconditionSelector<...>(D);
 *
 * // Inner solver: Typically limited to few iterations
 * //               using IterationNumberControl
 * auto S_inv_approx = inverse_operator(S_approx,...,S_approx_prec);
 *
 * // Construction of exact inverse of Schur complement
 * const auto S = schur_complement(A_inv,B,C,D);
 *
 * // Outer solver
 * const auto S_inv = inverse_operator(S,...,S_inv_approx);
 *
 * // Solve reduced block system
 * auto rhs = condense_schur_rhs (A_inv,C,f,g);
 *
 * // Solve for y
 * y = S_inv * rhs;
 * x = postprocess_schur_solution (A_inv,B,y,f);
 * @endcode
 * Note that due to the construction of @c S_inv_approx and subsequently @c
 * S_inv, there are a pair of nested iterative solvers which could
 * collectively consume a lot of resources. Therefore care should be taken in
 * the choices leading to the construction of the iterative inverse_operators.
 * One might consider the use of a IterationNumberControl (or a similar
 * mechanism) to limit the number of inner solver iterations. This controls
 * the accuracy of the approximate inverse operation $ \tilde{S}^{-1} $ which
 * acts only as the preconditioner for $ S^{-1} $. Furthermore, the
 * preconditioner to $ \tilde{S}^{-1} $, which in this example is $
 * \text{prec}(D) $, should ideally be computationally inexpensive.
 *
 * However, if an iterative solver based on IterationNumberControl is used as
 * a preconditioner then the preconditioning operation is not a linear
 * operation. Here a flexible solver like SolverFGMRES (flexible GMRES) is
 * best employed as an outer solver in order to deal with the variable
 * behavior of the preconditioner. Otherwise the iterative solver can
 * stagnate somewhere near the tolerance of the preconditioner or generally
 * behave erratically. Alternatively, using a ReductionControl would ensure
 * that the preconditioner always solves to the same tolerance, thereby
 * rendering its behavior constant.
 *
 * Further examples of this functionality can be found in the test-suite, such
 * as <code>tests/lac/schur_complement_01.cc</code>. The solution of a
 * multi-component problem (namely step-22) using the schur_complement can be
 * found in <code>tests/lac/schur_complement_03.cc</code>.
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 *
 * @ingroup LAOperators
 */
template <typename Range_1,
          typename Domain_1,
          typename Range_2,
          typename Domain_2,
          typename Payload>
LinearOperator<Range_2, Domain_2, Payload>
schur_complement(const LinearOperator<Domain_1, Range_1, Payload> &A_inv,
                 const LinearOperator<Range_1, Domain_2, Payload> &B,
                 const LinearOperator<Range_2, Domain_1, Payload> &C,
                 const LinearOperator<Range_2, Domain_2, Payload> &D)
{
  // We return the result of the compound LinearOperator
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  if (D.is_null_operator == false)
    return D - C * A_inv * B;
  else
    return -1.0 * C * A_inv * B;
}

/** @} */


/**
 * @name Creation of PackagedOperation objects related to the Schur Complement
 */
/** @{ */

/**
 * @relatesalso PackagedOperation
 *
 * For the system of equations
 * @f{eqnarray*}{
 *   Ax + By &=& f \\
 *   Cx + Dy &=& g \quad ,
 * @f}
 * this operation performs the pre-processing (condensation) step on the RHS
 * subvector @p g so that the Schur complement can be used to solve this
 * system of equations. More specifically, it produces an object that
 * represents the condensed form of the subvector @p g, namely
 * @f[
 *   g' = g - C \: A^{-1} \: f
 * @f]
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 *
 * @ingroup LAOperators
 */
template <typename Range_1,
          typename Domain_1,
          typename Range_2,
          typename Payload>
PackagedOperation<Range_2>
condense_schur_rhs(const LinearOperator<Range_1, Domain_1, Payload> &A_inv,
                   const LinearOperator<Range_2, Domain_1, Payload> &C,
                   const Range_1                                    &f,
                   const Range_2                                    &g)
{
  // We return the result of the compound PackagedOperation
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  return g - C * A_inv * f;
}

/**
 * @relatesalso PackagedOperation
 *
 * For the system of equations
 * @f{eqnarray*}{
 *   Ax + By &=& f \\
 *   Cx + Dy &=& g \quad ,
 * @f}
 * this operation performs the post-processing step of the Schur complement to
 * solve for the second subvector @p x once subvector @p y is known, with the
 * result that
 * @f[
 *   x =  A^{-1}(f - By)
 * @f]
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 *
 * @ingroup LAOperators
 */
template <typename Range_1,
          typename Domain_1,
          typename Domain_2,
          typename Payload>
PackagedOperation<Domain_1>
postprocess_schur_solution(
  const LinearOperator<Range_1, Domain_1, Payload> &A_inv,
  const LinearOperator<Range_1, Domain_2, Payload> &B,
  const Domain_2                                   &y,
  const Range_1                                    &f)
{
  // We return the result of the compound PackagedOperation
  // directly, so as to ensure that the underlying Payload
  // definition aligns with the operations expressed here.
  // All of the memory allocations etc. are taken care of
  // internally.
  return A_inv * (f - B * y);
}

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
