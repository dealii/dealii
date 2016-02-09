// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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

#ifndef dealii__mg_matrix_h
#define dealii__mg_matrix_h

#include <deal.II/lac/vector.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/

namespace mg
{
  /**
   * Multilevel matrix. This matrix stores an MGLevelObject of
   * PointerMatrixBase objects. It implements the interface defined in
   * MGMatrixBase, so that it can be used as a matrix in Multigrid.
   *
   * @author Guido Kanschat
   * @date 2002, 2010
   */
  template <typename VectorType = Vector<double> >
  class Matrix
    : public MGMatrixBase<VectorType>
  {
  public:
    /**
     * Default constructor for an empty object.
     */
    Matrix();

    /**
     * Constructor setting up pointers to the matrices in <tt>M</tt> by
     * calling initialize().
     */
    template <typename MatrixType>
    Matrix(const MGLevelObject<MatrixType> &M);

    /**
     * Initialize the object such that the level multiplication uses the
     * matrices in <tt>M</tt>
     */
    template <typename MatrixType>
    void
    initialize(const MGLevelObject<MatrixType> &M);

    /**
     * Access matrix on a level.
     */
    const PointerMatrixBase<VectorType> &operator[] (unsigned int level) const;

    virtual void vmult (const unsigned int level, VectorType &dst, const VectorType &src) const;
    virtual void vmult_add (const unsigned int level, VectorType &dst, const VectorType &src) const;
    virtual void Tvmult (const unsigned int level, VectorType &dst, const VectorType &src) const;
    virtual void Tvmult_add (const unsigned int level, VectorType &dst, const VectorType &src) const;

    /**
     * Memory used by this object.
     */
    std::size_t memory_consumption () const;
  private:
    MGLevelObject<std_cxx11::shared_ptr<PointerMatrixBase<VectorType> > > matrices;
  };

}


/**
 * Multilevel matrix selecting from block matrices. This class implements the
 * interface defined by MGMatrixBase.  The template parameter @p MatrixType
 * should be a block matrix class like BlockSparseMatrix or @p
 * BlockSparseMatrixEZ. Then, this class stores a pointer to a MGLevelObject
 * of this matrix class. In each @p vmult, the block selected on
 * initialization will be multiplied with the vector provided.
 *
 * @author Guido Kanschat, 2002
 */
template <typename MatrixType, typename number>
class MGMatrixSelect : public MGMatrixBase<Vector<number> >
{
public:
  /**
   * Constructor. @p row and @p col are the coordinate of the selected block.
   * The other argument is handed over to the @p SmartPointer constructor.
   */
  MGMatrixSelect (const unsigned int         row    = 0,
                  const unsigned int         col    = 0,
                  MGLevelObject<MatrixType> *matrix = 0);

  /**
   * Set the matrix object to be used. The matrix object must exist longer as
   * the @p MGMatrix object, since only a pointer is stored.
   */
  void set_matrix (MGLevelObject<MatrixType> *M);

  /**
   * Select the block for multiplication.
   */
  void select_block (const unsigned int row,
                     const unsigned int col);

  /**
   * Matrix-vector-multiplication on a certain level.
   */
  virtual void vmult (const unsigned int level,
                      Vector<number> &dst,
                      const Vector<number> &src) const;

  /**
   * Adding matrix-vector-multiplication on a certain level.
   */
  virtual void vmult_add (const unsigned int level,
                          Vector<number> &dst,
                          const Vector<number> &src) const;

  /**
   * Transpose matrix-vector-multiplication on a certain level.
   */
  virtual void Tvmult (const unsigned int level,
                       Vector<number> &dst,
                       const Vector<number> &src) const;

  /**
   * Adding transpose matrix-vector-multiplication on a certain level.
   */
  virtual void Tvmult_add (const unsigned int level,
                           Vector<number> &dst,
                           const Vector<number> &src) const;

private:
  /**
   * Pointer to the matrix objects on each level.
   */
  SmartPointer<MGLevelObject<MatrixType>,MGMatrixSelect<MatrixType,number> > matrix;
  /**
   * Row coordinate of selected block.
   */
  unsigned int row;
  /**
   * Column coordinate of selected block.
   */
  unsigned int col;

};

/*@}*/

/*----------------------------------------------------------------------*/

namespace mg
{
  template <typename VectorType>
  template <typename MatrixType>
  inline
  void
  Matrix<VectorType>::initialize (const MGLevelObject<MatrixType> &p)
  {
    matrices.resize(p.min_level(), p.max_level());
    for (unsigned int level=p.min_level(); level <= p.max_level(); ++level)
      matrices[level] = std_cxx11::shared_ptr<PointerMatrixBase<VectorType> >
                        (new_pointer_matrix_base(p[level], VectorType()));
  }


  template <typename VectorType>
  template <typename MatrixType>
  inline
  Matrix<VectorType>::Matrix (const MGLevelObject<MatrixType> &p)
  {
    initialize(p);
  }


  template <typename VectorType>
  inline
  Matrix<VectorType>::Matrix ()
  {}


  template <typename VectorType>
  inline
  const PointerMatrixBase<VectorType> &
  Matrix<VectorType>::operator[] (unsigned int level) const
  {
    return *matrices[level];
  }


  template <typename VectorType>
  void
  Matrix<VectorType>::vmult (const unsigned int level,
                             VectorType         &dst,
                             const VectorType   &src) const
  {
    matrices[level]->vmult(dst, src);
  }


  template <typename VectorType>
  void
  Matrix<VectorType>::vmult_add (const unsigned int level,
                                 VectorType         &dst,
                                 const VectorType   &src) const
  {
    matrices[level]->vmult_add(dst, src);
  }


  template <typename VectorType>
  void
  Matrix<VectorType>::Tvmult (const unsigned int level,
                              VectorType         &dst,
                              const VectorType   &src) const
  {
    matrices[level]->Tvmult(dst, src);
  }


  template <typename VectorType>
  void
  Matrix<VectorType>::Tvmult_add (const unsigned int level,
                                  VectorType         &dst,
                                  const VectorType   &src) const
  {
    matrices[level]->Tvmult_add(dst, src);
  }


  template <typename VectorType>
  inline
  std::size_t
  Matrix<VectorType>::memory_consumption () const
  {
    return sizeof(*this) + matrices->memory_consumption();
  }
}


/*----------------------------------------------------------------------*/

template <typename MatrixType, typename number>
MGMatrixSelect<MatrixType, number>::
MGMatrixSelect (const unsigned int         row,
                const unsigned int         col,
                MGLevelObject<MatrixType> *p)
  :
  matrix (p, typeid(*this).name()),
  row(row),
  col(col)
{}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::set_matrix (MGLevelObject<MatrixType> *p)
{
  matrix = p;
}


template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::
select_block (const unsigned int brow,
              const unsigned int bcol)
{
  row = brow;
  col = bcol;
}


template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::
vmult  (const unsigned int    level,
        Vector<number>       &dst,
        const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).vmult(dst, src);
}


template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::
vmult_add  (const unsigned int    level,
            Vector<number>       &dst,
            const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).vmult_add(dst, src);
}


template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::
Tvmult  (const unsigned int    level,
         Vector<number>       &dst,
         const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).Tvmult(dst, src);
}


template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::
Tvmult_add  (const unsigned int    level,
             Vector<number>       &dst,
             const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).Tvmult_add(dst, src);
}

DEAL_II_NAMESPACE_CLOSE

#endif
