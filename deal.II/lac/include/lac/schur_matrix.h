//-------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------
#ifndef __deal2__schur_matrix_h
#define __deal2__schur_matrix_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <base/logstream.h>
#include <lac/vector_memory.h>
#include <vector>

template <typename> class BlockVector;


/**
 * Schur complement of a block matrix.
 *
 * Given a non-singular matrix @p{A} (often positive definite) and a
 * positive semi-definite matrix @p{C} as well as matrices @p{B} and
 * @p{Dt} of full rank, this class implements a new matrix, the Schur
 * complement a the system of equations of the structure
 *
 * @begin{verbatim}
 * /        \  /   \     /   \
 * |  A  Dt |  | u |  -  | f |
 * | -B  C  |  | p |  -  | g |
 * \        /  \   /     \   /
 * @end{verbatim}
 *
 * Multiplication with the Schur matrix @p{S} is the operation
 * @begin{verbatim}
 * S p = C p + B A-inverse Dt-transpose p,
 * @end{verbatim}
 * which is an operation within the space for @p{p}.
 *
 * The data handed to the Schur matrix are as follows:
 *
 * @p{A}: the inverse of @p{A} is stored, instead of @p{A}. This
 * allows the application to use the most efficient form of inversion,
 * iterative or direct.
 *
 * @p{B}, @p{C}: these matrices are stored "as is".
 *
 * @p{Dt}: the computation of the Schur complement involves the
 * function @p{Tvmult} of the matrix @p{Dt}, not @p{vmult}! This way,
 * it is sufficient to build only one matrix @p{B} for the symmetric
 * Schur complement and use it twice.
 *
 * All matrices involved are of arbitrary type and vectors are
 * @ref{BlockVector}s. This way, @p{SchurMatrix} can be coupled with
 * any matrix classes providing @p{vmult} and @p{Tvmult} and can be
 * even nested. Since @ref{SmartPointer}s of matrices are stored, the
 * matrix blocks should be derived from @ref{Subscriptor}.
 *
 * Since the Schur complement of a matrix corresponds to a Gaussian
 * block elimination, the right hand side of the condensed system must
 * be preprocessed. Furthermore, the eliminated variable must be
 * reconstructed after solving.
 *
 * @begin{verbatim}
 *   g = g + B A-inverse f
 *   u = A-inverse (f - D-transpose p)
 * @end{verbatim}
 *
 * Applying these transformations, the solution of the system above by a
 * @p{SchurMatrix} @p{schur} is coded as follows:
 *
 * @begin{verbatim}
 *   schur.prepare_rhs (g, f);
 *   solver.solve (schur, p, g, precondition);
 *   schur.postprocess (u, p);
 * @end{verbatim}
 *
 * @author Guido Kanschat, 2000, 2001
 */
template <class MA_inverse, class MB, class MDt, class MC>
class SchurMatrix :
  public Subscriptor
{
  public:
  SchurMatrix(const MA_inverse& Ainv,
	      const MB& B,
	      const MDt& Dt,
	      const MC& C,
	      VectorMemory<BlockVector<double> >& mem);

				   /**
				    * Do block elimination of the
				    * right hand side. Given right
				    * hand sides for both components
				    * of the block system, this
				    * function provides the right hand
				    * side for the Schur complement.
				    *
				    * The result is stored in the
				    * first argument, which is also
				    * part of the input data. If it is
				    * necessary to conserve the data,
				    * @p{dst} must be copied before
				    * calling this function. This is
				    * reasonable, since in many cases,
				    * only the pre-processed right
				    * hand side is needed.
				    */
  void prepare_rhs (BlockVector<double>& dst,
		    const BlockVector<double>& src) const;

				   /**
				    * Multiplication with the Schur
				    * complement.
				    */
  void vmult (BlockVector<double>& dst,
	      const BlockVector<double>& src) const;

//  void Tmult(BlockVector<double>& dst, const BlockVector<double>& src) const;

				   /**
				    * Computation of the residual of
				    * the Schur complement.
				    */
  double residual (BlockVector<double>& dst,
		   const BlockVector<double>& src,
		   const BlockVector<double>& rhs) const;

				   /**
				    * Compute the eliminated variable
				    * from the solution of the Schur
				    * complement problem.
				    */
  void postprocess (BlockVector<double>& dst,
		    const BlockVector<double>& src,
		    const BlockVector<double>& rhs) const;
  private:
				   /**
				    * No copy constructor.
				    */
  SchurMatrix (const SchurMatrix<MA_inverse, MB, MDt, MC>&);
				   /**
				    * No assignment.
				    */
  SchurMatrix& operator = (const SchurMatrix<MA_inverse, MB, MDt, MC>&);

				   /**
				    * Pointer to inverse of upper left block.
				    */
  const SmartPointer<const MA_inverse> Ainv;
				   /**
				    * Pointer to lower left block.
				    */
  const SmartPointer<const MB> B;
				   /**
				    * Pointer to transpose of upper right block.
				    */
  const SmartPointer<const MDt> Dt;
				   /**
				    * Pointer to lower right block.
				    */
  const SmartPointer<const MC> C;
				   /**
				    * Auxiliary memory for vectors.
				    */
  VectorMemory<BlockVector<double> >& mem;
};

template <class MA_inverse, class MB, class MDt, class MC>
SchurMatrix<MA_inverse, MB, MDt, MC>
::SchurMatrix(const MA_inverse& Ainv,
	      const MB& B,
	      const MDt& Dt,
	      const MC& C,
	      VectorMemory<BlockVector<double> >& mem)
  : Ainv(&Ainv), B(&B), Dt(&Dt), C(&C), mem(mem)
{
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::vmult(BlockVector<double>& dst,
	const BlockVector<double>& src) const
{
  deallog.push("Schur");
  C->vmult(dst, src);
  BlockVector<double>* h1 = mem.alloc();
  h1->reinit(B->n_block_cols(), src.block(0).size());
  Dt->Tvmult(*h1,src);
  BlockVector<double>* h2 = mem.alloc();
  h2->reinit(*h1);
  Ainv->vmult(*h2, *h1);
  mem.free(h1);
  B->vmult_add(dst, *h2);
  mem.free(h2);
  deallog.pop();
}


template <class MA_inverse, class MB, class MDt, class MC>
double SchurMatrix<MA_inverse, MB, MDt, MC>
::residual(BlockVector<double>& dst,
	   const BlockVector<double>& src,
	   const BlockVector<double>& rhs) const
{
  vmult(dst, src);
  dst.scale(-1.);
  dst += rhs;
  return dst.l2_norm();
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::prepare_rhs(BlockVector<double>& dst,
	      const BlockVector<double>& src) const
{
  Assert (src.n_blocks() == B->n_block_cols(),
	  ExcDimensionMismatch(src.n_blocks(), B->n_block_cols()));
  Assert (dst.n_blocks() == B->n_block_rows(),
	  ExcDimensionMismatch(dst.n_blocks(), B->n_block_rows()));
  
  deallog.push("Schur-prepare");
  BlockVector<double>* h1 = mem.alloc();
  h1->reinit(B->n_block_cols(), src.block(0).size());
  Ainv->vmult(*h1, src);
  B->vmult_add(dst, *h1);
  mem.free(h1);
  deallog.pop();
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::postprocess(BlockVector<double>& dst,
	      const BlockVector<double>& src,
	      const BlockVector<double>& rhs) const
{
  Assert (dst.n_blocks() == B->n_block_cols(),
	  ExcDimensionMismatch(dst.n_blocks(), B->n_block_cols()));
  Assert (rhs.n_blocks() == B->n_block_cols(),
	  ExcDimensionMismatch(rhs.n_blocks(), B->n_block_cols()));
  Assert (src.n_blocks() == B->n_block_rows(),
	  ExcDimensionMismatch(src.n_blocks(), B->n_block_rows()));
  
  deallog.push("Schur-post");
  BlockVector<double>* h1 = mem.alloc();
  h1->reinit(B->n_block_cols(), src.block(0).size());
  Dt->Tvmult(*h1, src);
  h1->sadd(-1.,rhs);
  Ainv->vmult(dst,*h1);
  mem.free(h1);
  deallog.pop();
}


#endif
