// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2016 by the deal.II authors
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

#ifndef dealii__matrix_lib_h
#define dealii__matrix_lib_h

#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/solver_richardson.h>

DEAL_II_NAMESPACE_OPEN

template<typename number> class Vector;
template<typename number> class BlockVector;
template<typename number> class SparseMatrix;

/*! @addtogroup Matrix2
 *@{
 */


/**
 * A matrix that is the multiple of another matrix.
 *
 * Matrix-vector products of this matrix are composed of those of the original
 * matrix with the vector and then scaling of the result by a constant factor.
 *
 * @deprecated If deal.II was configured with C++11 support, use the
 * LinearOperator class instead, see the module on
 * @ref LAOperators "linear operators"
 * for further details.
 *
 * @author Guido Kanschat, 2007
 */
template<typename VectorType>
class ScaledMatrix : public Subscriptor
{
public:
  /**
   * Constructor leaving an uninitialized object.
   */
  ScaledMatrix ();
  /**
   * Constructor with initialization.
   */
  template <typename MatrixType>
  ScaledMatrix (const MatrixType &M,
                const double factor);

  /**
   * Destructor
   */
  ~ScaledMatrix ();

  /**
   * Initialize for use with a new matrix and factor.
   */
  template <typename MatrixType>
  void initialize (const MatrixType &M,
                   const double factor);

  /**
   * Reset the object to its original state.
   */
  void clear ();

  /**
   * Matrix-vector product.
   */
  void vmult (VectorType &w,
              const VectorType &v) const;

  /**
   * Transposed matrix-vector product.
   */
  void Tvmult (VectorType &w,
               const VectorType &v) const;

private:
  /**
   * The matrix.
   */
  PointerMatrixBase<VectorType> *m;

  /**
   * The scaling factor;
   */
  double factor;
};


/**
 * Mean value filter.  The vmult() functions of this matrix filter out mean
 * values of the vector.  If the vector is of type BlockVector, then an
 * additional parameter selects a single component for this operation.
 *
 * In mathematical terms, this class acts as if it was the matrix $I-\frac
 * 1n{\mathbf 1}_n{\mathbf 1}_n^T$ where ${\mathbf 1}_n$ is a vector of size
 * $n$ that has only ones as its entries. Thus, taking the dot product between
 * a vector $\mathbf v$ and $\frac 1n {\mathbf 1}_n$ yields the <i>mean
 * value</i> of the entries of ${\mathbf v}$. Consequently, $ \left[I-\frac
 * 1n{\mathbf 1}_n{\mathbf 1}_n^T\right] \mathbf v = \mathbf v - \left[\frac
 * 1n {\mathbf v} \cdot {\mathbf 1}_n\right]{\mathbf 1}_n$ subtracts from every
 * vector element the mean value of all elements.
 *
 * @author Guido Kanschat, 2002, 2003
 */
class MeanValueFilter : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor, optionally selecting a component.
   */
  MeanValueFilter(const size_type component = numbers::invalid_size_type);

  /**
   * Subtract mean value from @p v.
   */
  template <typename number>
  void filter (Vector<number> &v) const;

  /**
   * Subtract mean value from @p v.
   */
  template <typename number>
  void filter (BlockVector<number> &v) const;

  /**
   * Return the source vector with subtracted mean value.
   */
  template <typename number>
  void vmult (Vector<number>       &dst,
              const Vector<number> &src) const;

  /**
   * Add source vector with subtracted mean value to dest.
   */
  template <typename number>
  void vmult_add (Vector<number>       &dst,
                  const Vector<number> &src) const;

  /**
   * Return the source vector with subtracted mean value in selected
   * component.
   */
  template <typename number>
  void vmult (BlockVector<number>       &dst,
              const BlockVector<number> &src) const;

  /**
   * Add a source to dest, where the mean value in the selected component is
   * subtracted.
   */
  template <typename number>
  void vmult_add (BlockVector<number>       &dst,
                  const BlockVector<number> &src) const;


  /**
   * Not implemented.
   */
  template <typename VectorType>
  void Tvmult(VectorType &, const VectorType &) const;

  /**
   * Not implemented.
   */
  template <typename VectorType>
  void Tvmult_add(VectorType &, const VectorType &) const;

private:
  /**
   * Component for filtering block vectors.
   */
  const size_type component;
};



/**
 * Objects of this type represent the inverse of a matrix as computed
 * approximately by using the SolverRichardson iterative solver. In other
 * words, if you set up an object of the current type for a matrix $A$, then
 * calling the vmult() function with arguments $v,w$ amounts to setting
 * $w=A^{-1}v$ by solving the linear system $Aw=v$ using the Richardson solver
 * with a preconditioner that can be chosen. Similarly, this class allows to
 * also multiple with the transpose of the inverse (i.e., the inverse of the
 * transpose) using the function SolverRichardson::Tsolve().
 *
 * The functions vmult() and Tvmult() approximate the inverse iteratively
 * starting with the vector <tt>dst</tt>. Functions vmult_add() and
 * Tvmult_add() start the iteration with a zero vector. All of the matrix-
 * vector multiplication functions expect that the Richardson solver with the
 * given preconditioner actually converge. If the Richardson solver does not
 * converge within the specified number of iterations, the exception that will
 * result in the solver will simply be propagated to the caller of the member
 * function of the current class.
 *
 * @note A more powerful version of this class is provided by the
 * IterativeInverse class.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on
 * @ref Instantiations
 * in the manual).
 *
 * @deprecated If deal.II was configured with C++11 support, use the
 * LinearOperator class instead, see the module on
 * @ref LAOperators "linear operators"
 * for further details.
 *
 * @author Guido Kanschat, 2005
 */
template<typename VectorType>
class InverseMatrixRichardson : public Subscriptor
{
public:
  /**
   * Constructor, initializing the solver with a control and memory object.
   * The inverted matrix and the preconditioner are added in initialize().
   */
  InverseMatrixRichardson (SolverControl            &control,
                           VectorMemory<VectorType> &mem);
  /**
   * Since we use two pointers, we must implement a destructor.
   */
  ~InverseMatrixRichardson();

  /**
   * Initialization function. Provide a solver object, a matrix, and another
   * preconditioner for this.
   */
  template <typename MatrixType, typename PreconditionerType>
  void initialize (const MatrixType &,
                   const PreconditionerType &);

  /**
   * Access to the SolverControl object used by the solver.
   */
  SolverControl &control() const;
  /**
   * Execute solver.
   */
  void vmult (VectorType &, const VectorType &) const;

  /**
   * Execute solver.
   */
  void vmult_add (VectorType &, const VectorType &) const;

  /**
   * Execute transpose solver.
   */
  void Tvmult (VectorType &, const VectorType &) const;

  /**
   * Execute transpose solver.
   */
  void Tvmult_add (VectorType &, const VectorType &) const;

private:
  /**
   * A reference to the provided VectorMemory object.
   */
  VectorMemory<VectorType> &mem;

  /**
   * The solver object.
   */
  mutable SolverRichardson<VectorType> solver;

  /**
   * The matrix in use.
   */
  PointerMatrixBase<VectorType> *matrix;

  /**
   * The preconditioner to use.
   */
  PointerMatrixBase<VectorType> *precondition;
};




/*@}*/
//---------------------------------------------------------------------------


template<typename VectorType>
inline
ScaledMatrix<VectorType>::ScaledMatrix()
  :
  m(0),
  factor (0)
{}



template<typename VectorType>
template<typename MatrixType>
inline
ScaledMatrix<VectorType>::ScaledMatrix(const MatrixType &mat, const double factor)
  :
  m(new_pointer_matrix_base(mat, VectorType())),
  factor(factor)
{}



template<typename VectorType>
template<typename MatrixType>
inline
void
ScaledMatrix<VectorType>::initialize(const MatrixType &mat, const double f)
{
  if (m) delete m;
  m = new_pointer_matrix_base(mat, VectorType());
  factor = f;
}



template<typename VectorType>
inline
void
ScaledMatrix<VectorType>::clear()
{
  if (m) delete m;
  m = 0;
}



template<typename VectorType>
inline
ScaledMatrix<VectorType>::~ScaledMatrix()
{
  clear ();
}


template<typename VectorType>
inline
void
ScaledMatrix<VectorType>::vmult (VectorType &w, const VectorType &v) const
{
  m->vmult(w, v);
  w *= factor;
}


template<typename VectorType>
inline
void
ScaledMatrix<VectorType>::Tvmult (VectorType &w, const VectorType &v) const
{
  m->Tvmult(w, v);
  w *= factor;
}


//---------------------------------------------------------------------------

template <typename VectorType>
inline void
MeanValueFilter::Tvmult(VectorType &, const VectorType &) const
{
  Assert(false, ExcNotImplemented());
}


template <typename VectorType>
inline void
MeanValueFilter::Tvmult_add(VectorType &, const VectorType &) const
{
  Assert(false, ExcNotImplemented());
}

//-----------------------------------------------------------------------//

template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
inline void
InverseMatrixRichardson<VectorType>::initialize (const MatrixType &m,
                                                 const PreconditionerType &p)
{
  if (matrix != 0)
    delete matrix;
  matrix = new PointerMatrix<MatrixType, VectorType>(&m);
  if (precondition != 0)
    delete precondition;
  precondition = new PointerMatrix<PreconditionerType, VectorType>(&p);
}


DEAL_II_NAMESPACE_CLOSE

#endif
