// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__iterative_inverse_h
#define __deal2__iterative_inverse_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/pointer_matrix.h>

DEAL_II_NAMESPACE_OPEN


/**
 * Implementation of the inverse of a matrix, using an iterative
 * method.
 *
 * The function vmult() of this class starts an iterative solver in
 * order to approximate the action of the inverse matrix.
 *
 * Krylov space methods like SolverCG or SolverBicgstab
 * become inefficient if solution down to machine accuracy is
 * needed. This is due to the fact, that round-off errors spoil the
 * orthogonality of the vector sequences. Therefore, a nested
 * iteration of two methods is proposed: The outer method is
 * SolverRichardson, since it is robust with respect to round-of
 * errors. The inner loop is an appropriate Krylov space method, since
 * it is fast.
 *
 * @code
 * // Declare related objects
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * GrowingVectorMemory<Vector<double> > mem;
 *
 * ReductionControl inner_control (10, 1.e-30, 1.e-2);
 * PreconditionSSOR <SparseMatrix<double> > inner_precondition;
 * inner_precondition.initialize (A, 1.2);
 *
 * IterativeInverse<Vector<double> > precondition;
 * precondition.initialize (A, inner_precondition);
 * precondition.solver.select("cg");
 * precondition.solver.set_control(inner_control);
 *
 * SolverControl outer_control(100, 1.e-16);
 * SolverRichardson<Vector<double> > outer_iteration;
 *
 * outer_iteration.solve (A, x, b, precondition);
 * @endcode
 *
 * Each time we call the inner loop, reduction of the residual by a
 * factor <tt>1.e-2</tt> is attempted. Since the right hand side vector of
 * the inner iteration is the residual of the outer loop, the relative
 * errors are far from machine accuracy, even if the errors of the
 * outer loop are in the range of machine accuracy.
 *
 * @ingroup Matrix2
 * @author Guido Kanschat
 * @date 2010
 */
template <class VECTOR>
class IterativeInverse : public Subscriptor
{
public:
  /**
   * Initialization
   * function. Provide a matrix and
   * preconditioner for the solve in
   * vmult().
   */
  template <class MATRIX, class PRECONDITION>
  void initialize (const MATRIX &, const PRECONDITION &);

  /**
   * Delete the pointers to matrix
   * and preconditioner.
   */
  void clear();

  /**
   * Solve for right hand side <tt>src</tt>.
   */
  void vmult (VECTOR &dst, const VECTOR &src) const;

  /**
   * Solve for right hand side <tt>src</tt>, but allow for the fact
   * that the vectors given to this function have different type from
   * the vectors used by the inner solver.
   */
  template <class VECTOR2>
  void vmult (VECTOR2 &dst, const VECTOR2 &src) const;

  /**
   * The solver, which allows
   * selection of the actual solver
   * as well as adjuxtment of
   * parameters.
   */
  SolverSelector<VECTOR> solver;

private:
  /**
   * The matrix in use.
   */
  std_cxx11::shared_ptr<PointerMatrixBase<VECTOR> > matrix;

  /**
   * The preconditioner to use.
   */
  std_cxx11::shared_ptr<PointerMatrixBase<VECTOR> > preconditioner;
};


template <class VECTOR>
template <class MATRIX, class PRECONDITION>
inline
void
IterativeInverse<VECTOR>::initialize(const MATRIX &m, const PRECONDITION &p)
{
  // dummy variable
  VECTOR *v = 0;
  matrix = std_cxx11::shared_ptr<PointerMatrixBase<VECTOR> > (new_pointer_matrix_base(m, *v));
  preconditioner = std_cxx11::shared_ptr<PointerMatrixBase<VECTOR> > (new_pointer_matrix_base(p, *v));
}


template <class VECTOR>
inline
void
IterativeInverse<VECTOR>::clear()
{
  matrix = 0;
  preconditioner = 0;
}


template <class VECTOR>
inline void
IterativeInverse<VECTOR>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert(matrix.get() != 0, ExcNotInitialized());
  Assert(preconditioner.get() != 0, ExcNotInitialized());
  dst = 0.;
  solver.solve(*matrix, dst, src, *preconditioner);
}


template <class VECTOR>
template <class VECTOR2>
inline void
IterativeInverse<VECTOR>::vmult (VECTOR2 &dst, const VECTOR2 &src) const
{
  Assert(matrix.get() != 0, ExcNotInitialized());
  Assert(preconditioner.get() != 0, ExcNotInitialized());
  GrowingVectorMemory<VECTOR> mem;
  typename VectorMemory<VECTOR>::Pointer sol(mem);
  typename VectorMemory<VECTOR>::Pointer rhs(mem);
  sol->reinit(dst);
  *rhs = src;
  solver.solve(*matrix, *sol, *rhs, *preconditioner);
  dst = *sol;
}



DEAL_II_NAMESPACE_CLOSE

#endif


