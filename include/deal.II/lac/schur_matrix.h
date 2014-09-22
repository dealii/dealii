// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__schur_matrix_h
#define __deal2__schur_matrix_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/block_vector.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Matrix2
 *@{
 */

/**
 * Schur complement of a block matrix.
 *
 * Given a non-singular matrix @p A (often positive definite) and a
 * positive semi-definite matrix @p C as well as matrices @p B and
 * @p Dt of full rank, this class implements a new matrix, the Schur
 * complement a the system of equations of the structure
 *
 * @verbatim
 * /        \  /   \     /   \
 * |  A  Dt |  | u |  -  | f |
 * | -B  C  |  | p |  -  | g |
 * \        /  \   /     \   /
 * @endverbatim
 *
 * Multiplication with the Schur matrix @p S is the operation
 * @verbatim
 * S p = C p + B A-inverse Dt-transpose p,
 * @endverbatim
 * which is an operation within the space for @p p.
 *
 * The data handed to the Schur matrix are as follows:
 *
 * @p A: the inverse of @p A is stored, instead of @p A. This
 * allows the application to use the most efficient form of inversion,
 * iterative or direct.
 *
 * @p B, @p C: these matrices are stored "as is".
 *
 * @p Dt: the computation of the Schur complement involves the
 * function @p Tvmult of the matrix @p Dt, not @p vmult! This way,
 * it is sufficient to build only one matrix @p B for the symmetric
 * Schur complement and use it twice.
 *
 * All matrices involved are of arbitrary type and vectors are
 * BlockVectors. This way, @p SchurMatrix can be coupled with
 * any matrix classes providing @p vmult and @p Tvmult and can be
 * even nested. Since SmartPointers of matrices are stored, the
 * matrix blocks should be derived from Subscriptor.
 *
 * Since the Schur complement of a matrix corresponds to a Gaussian
 * block elimination, the right hand side of the condensed system must
 * be preprocessed. Furthermore, the eliminated variable must be
 * reconstructed after solving.
 *
 * @verbatim
 *   g = g + B A-inverse f
 *   u = A-inverse (f - D-transpose p)
 * @endverbatim
 *
 * Applying these transformations, the solution of the system above by a
 * @p SchurMatrix @p schur is coded as follows:
 *
 * @code
 *   schur.prepare_rhs (g, f);
 *   solver.solve (schur, p, g, precondition);
 *   schur.postprocess (u, p);
 * @endcode
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat, 2000, 2001, 2002
 */
template <class MA_inverse, class MB, class MDt, class MC>
class SchurMatrix : public Subscriptor
{
public:

  /**
   * Constructor. This constructor
   * receives all the matrices
   * needed. Furthermore, it gets a
   * reference to a memory structure
   * for obtaining block vectors.
   *
   * Optionally, the length of the
   * @p u-vector can be provided.
   *
   * For the meaning of the matrices
   * see the class documentation.
   */
  SchurMatrix(const MA_inverse &Ainv,
              const MB &B,
              const MDt &Dt,
              const MC &C,
              VectorMemory<BlockVector<double> > &mem,
              const std::vector<unsigned int> &signature = std::vector<unsigned int>(0));

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
   * @p dst must be copied before
   * calling this function. This is
   * reasonable, since in many cases,
   * only the pre-processed right
   * hand side is needed.
   */
  void prepare_rhs (BlockVector<double> &dst,
                    const BlockVector<double> &src) const;

  /**
   * Multiplication with the Schur
   * complement.
   */
  void vmult (BlockVector<double> &dst,
              const BlockVector<double> &src) const;

//  void Tmult(BlockVector<double>& dst, const BlockVector<double>& src) const;

  /**
   * Computation of the residual of
   * the Schur complement.
   */
  double residual (BlockVector<double> &dst,
                   const BlockVector<double> &src,
                   const BlockVector<double> &rhs) const;

  /**
   * Compute the eliminated variable
   * from the solution of the Schur
   * complement problem.
   */
  void postprocess (BlockVector<double> &dst,
                    const BlockVector<double> &src,
                    const BlockVector<double> &rhs) const;

  /**
   * Select debugging information for
   * log-file.  Debug level 1 is
   * defined and writes the norm of
   * every vector before and after
   * each operation. Debug level 0
   * turns off debugging information.
   */
  void debug_level(unsigned int l);
private:
  /**
   * No copy constructor.
   */
  SchurMatrix (const SchurMatrix<MA_inverse, MB, MDt, MC> &);
  /**
   * No assignment.
   */
  SchurMatrix &operator = (const SchurMatrix<MA_inverse, MB, MDt, MC> &);

  /**
   * Pointer to inverse of upper left block.
   */
  const SmartPointer<const MA_inverse,SchurMatrix<MA_inverse,MB,MDt,MC> > Ainv;
  /**
   * Pointer to lower left block.
   */
  const SmartPointer<const MB,SchurMatrix<MA_inverse,MB,MDt,MC> > B;
  /**
   * Pointer to transpose of upper right block.
   */
  const SmartPointer<const MDt,SchurMatrix<MA_inverse,MB,MDt,MC> > Dt;
  /**
   * Pointer to lower right block.
   */
  const SmartPointer<const MC,SchurMatrix<MA_inverse,MB,MDt,MC> > C;
  /**
   * Auxiliary memory for vectors.
   */
  VectorMemory<BlockVector<double> > &mem;

  /**
   * Optional signature of the @p u-vector.
   */
  std::vector<types::global_dof_index> signature;

  /**
   * Switch for debugging information.
   */
  unsigned int debug;
};

/*@}*/
//---------------------------------------------------------------------------

template <class MA_inverse, class MB, class MDt, class MC>
SchurMatrix<MA_inverse, MB, MDt, MC>
::SchurMatrix(const MA_inverse &Ainv,
              const MB &B,
              const MDt &Dt,
              const MC &C,
              VectorMemory<BlockVector<double> > &mem,
              const std::vector<unsigned int> &signature)
  : Ainv(&Ainv), B(&B), Dt(&Dt), C(&C),
    mem(mem),
    signature(signature),
    debug(0)
{
}


template <class MA_inverse, class MB, class MDt, class MC>
void
SchurMatrix<MA_inverse, MB, MDt, MC>
::debug_level(unsigned int l)
{
  debug = l;
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::vmult(BlockVector<double> &dst,
        const BlockVector<double> &src) const
{
  deallog.push("Schur");
  if (debug > 0)
    deallog << "src:" << src.l2_norm() << std::endl;

  C->vmult(dst, src);
  if (debug > 0)
    deallog << "C:" << dst.l2_norm() << std::endl;

  BlockVector<double> *h1 = mem.alloc();
  if (signature.size()>0)
    h1->reinit(signature);
  else
    h1->reinit(B->n_block_cols(), src.block(0).size());
  Dt->Tvmult(*h1,src);
  if (debug > 0)
    deallog << "Dt:" << h1->l2_norm() << std::endl;

  BlockVector<double> *h2 = mem.alloc();
  h2->reinit(*h1);
  Ainv->vmult(*h2, *h1);
  if (debug > 0)
    deallog << "Ainverse:" << h2->l2_norm() << std::endl;

  mem.free(h1);
  B->vmult_add(dst, *h2);
  if (debug > 0)
    deallog << "dst:" << dst.l2_norm() << std::endl;

  mem.free(h2);
  deallog.pop();
}


template <class MA_inverse, class MB, class MDt, class MC>
double SchurMatrix<MA_inverse, MB, MDt, MC>
::residual(BlockVector<double> &dst,
           const BlockVector<double> &src,
           const BlockVector<double> &rhs) const
{
  vmult(dst, src);
  dst.scale(-1.);
  dst += rhs;
  return dst.l2_norm();
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::prepare_rhs(BlockVector<double> &dst,
              const BlockVector<double> &src) const
{
  Assert (src.n_blocks() == B->n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), B->n_block_cols()));
  Assert (dst.n_blocks() == B->n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), B->n_block_rows()));

  deallog.push("Schur-prepare");
  if (debug > 0)
    deallog << "src:" << src.l2_norm() << std::endl;
  BlockVector<double> *h1 = mem.alloc();
  if (signature.size()>0)
    h1->reinit(signature);
  else
    h1->reinit(B->n_block_cols(), src.block(0).size());
  Ainv->vmult(*h1, src);
  if (debug > 0)
    deallog << "Ainverse:" << h1->l2_norm() << std::endl;
  B->vmult_add(dst, *h1);
  if (debug > 0)
    deallog << "dst:" << dst.l2_norm() << std::endl;
  mem.free(h1);
  deallog.pop();
}


template <class MA_inverse, class MB, class MDt, class MC>
void SchurMatrix<MA_inverse, MB, MDt, MC>
::postprocess(BlockVector<double> &dst,
              const BlockVector<double> &src,
              const BlockVector<double> &rhs) const
{
  Assert (dst.n_blocks() == B->n_block_cols(),
          ExcDimensionMismatch(dst.n_blocks(), B->n_block_cols()));
  Assert (rhs.n_blocks() == B->n_block_cols(),
          ExcDimensionMismatch(rhs.n_blocks(), B->n_block_cols()));
  Assert (src.n_blocks() == B->n_block_rows(),
          ExcDimensionMismatch(src.n_blocks(), B->n_block_rows()));

  deallog.push("Schur-post");
  if (debug > 0)
    deallog << "src:" << src.l2_norm() << std::endl;
  BlockVector<double> *h1 = mem.alloc();
  if (signature.size()>0)
    h1->reinit(signature);
  else
    h1->reinit(B->n_block_cols(), src.block(0).size());
  Dt->Tvmult(*h1, src);
  if (debug > 0)
    deallog << "Dt:" << h1->l2_norm() << std::endl;
  h1->sadd(-1.,rhs);
  Ainv->vmult(dst,*h1);
  if (debug > 0)
    deallog << "dst:" << dst.l2_norm() << std::endl;
  mem.free(h1);
  deallog.pop();
}


DEAL_II_NAMESPACE_CLOSE

#endif
