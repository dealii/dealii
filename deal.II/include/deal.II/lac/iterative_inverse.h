//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__iterative_inverse_h
#define __deal2__iterative_inverse_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
//#include <deal.II/base/template_constraints.h>
//#include <deal.II/lac/tridiagonal_matrix.h>
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
 * become inefficient if soution down to machine accuracy is
 * needed. This is due to the fact, that round-off errors spoil the
 * orthogonality of the vector sequences. Therefore, a nested
 * iteration of two methods is proposed: The outer method is
 * SolverRichardson, since it is robust with respect to round-of
 * errors. The inner loop is an appropriate Krylov space method, since
 * it is fast.
 *
 * @code
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * GrowingVectorMemory<Vector<double> > mem;

 * ReductionControl inner_control (10, 1.e-30, 1.e-2)
 * PreconditionSSOR <SparseMatrix<double> > inner_precondition;
 * inner_precondition.initialize (A, 1.2);
 *
 * IterativeInverse<Vector<double> > precondition;
 * precondition.initialize (A, inner_precondition);
 * precondition.solver.select("cg");
 * precondition.solver.control = inner_control;
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
				      * preconditioner for th solve in
				      * vmult().
				      */
    template <class MATRIX, class PRECONDITION>
    void initialize (const MATRIX&, const PRECONDITION&);

				     /**
				      * Delete the pointers to matrix
				      * and preconditioner.
				      */
    void clear();

				     /**
				      * Solve for right hand side <tt>src</tt>.
				      */
    void vmult (VECTOR& dst, const VECTOR& src) const;

				     /**
				      * Solve the transposed system
				      * for right hand side
				      * <tt>src</tt>.
				      */
//    void Tvmult (VECTOR& dst, const VECTOR& src) const;

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
    std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > matrix;

				     /**
				      * The preconditioner to use.
				      */
    std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > preconditioner;
				     /**
				      * The transpose of the matrix in use.
				      */
    std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > transpose_matrix;

				     /**
				      * The transpose of the preconditioner to use.
				      */
    std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > transpose_preconditioner;
};


template <class VECTOR>
template <class MATRIX, class PRECONDITION>
inline
void
IterativeInverse<VECTOR>::initialize(const MATRIX& m, const PRECONDITION& p)
{
				   // dummy variable
  VECTOR* v = 0;
  matrix = std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > (new_pointer_matrix_base(m, *v));
  preconditioner = std_cxx1x::shared_ptr<PointerMatrixBase<VECTOR> > (new_pointer_matrix_base(p, *v));
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
IterativeInverse<VECTOR>::vmult (VECTOR& dst, const VECTOR& src) const
{
  Assert(matrix.get() != 0, ExcNotInitialized());
  Assert(preconditioner.get() != 0, ExcNotInitialized());
  solver.solve(*matrix, dst, src, *preconditioner);
}


// template <class VECTOR>
// inline void
// IterativeInverse<VECTOR>::Tvmult (VECTOR& dst, const VECTOR& src) const
// {
//   TransposeMatrix<
//   solver.solve(*matrix, dst, src, *preconditioner);
// }


DEAL_II_NAMESPACE_CLOSE

#endif


