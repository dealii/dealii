// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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

#ifndef __deal2__pointer_matrix_h
#define __deal2__pointer_matrix_h

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template<class VECTOR> class VectorMemory;

class IdentityMatrix;
template <typename number> class FullMatrix;
template <typename number> class LAPACKFullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class BlockSparseMatrix;
template <typename number> class SparseMatrixEZ;
template <typename number> class BlockSparseMatrixEZ;
template <typename number> class BlockMatrixArray;
template <typename number> class TridiagonalMatrix;


/*! @addtogroup Matrix2
 *@{
 */

/**
 * Abstract class for use in iterations.  This class provides the
 * interface required by LAC solver classes. It allows to use
 * different concrete matrix classes in the same context, as long as
 * they apply to the same vector class.
 *
 * @author Guido Kanschat, 2000, 2001, 2002
 */
template<class VECTOR>
class PointerMatrixBase : public Subscriptor
{
public:
  /**
   * Value type of this
   * matrix. since the matrix
   * itself is unknown, we take the
   * value type of the
   * vector. Therefore, matrix
   * entries must be convertible to
   * vector entries.
   *
   * This was defined to make this
   * matrix a possible template
   * argument to
   * BlockMatrixArray.
   */
  typedef typename VECTOR::value_type value_type;

  /**
   * Virtual destructor.  Does
   * nothing except making sure that
   * the destructor of any derived
   * class is called whenever a pointer-to-base-class object is destroyed.
   */
  virtual ~PointerMatrixBase ();

  /**
   * Reset the object to its original state.
   */
  virtual void clear () = 0;

  /**
   * Matrix-vector product.
   */
  virtual void vmult (VECTOR &dst,
                      const VECTOR &src) const = 0;

  /**
   * Transposed matrix-vector product.
   */
  virtual void Tvmult (VECTOR &dst,
                       const VECTOR &src) const = 0;

  /**
   * Matrix-vector product, adding to
   * <tt>dst</tt>.
   */
  virtual void vmult_add (VECTOR &dst,
                          const VECTOR &src) const = 0;

  /**
   * Transposed matrix-vector product,
   * adding to <tt>dst</tt>.
   */
  virtual void Tvmult_add (VECTOR &dst,
                           const VECTOR &src) const = 0;
};


/**
 * A pointer to be used as a matrix.  This class stores a pointer to a
 * matrix and can be used as a matrix itself in iterative methods.
 *
 * The main purpose for the existence of this class is its base class,
 * which only has a vector as template argument. Therefore, this
 * interface provides an abstract base class for matrices.
 *
 * @author Guido Kanschat 2000, 2001, 2002
 */
template<class MATRIX, class VECTOR>
class PointerMatrix : public PointerMatrixBase<VECTOR>
{
public:
  /**
   * Constructor.  The pointer in the
   * argument is stored in this
   * class. As usual, the lifetime of
   * <tt>*M</tt> must be longer than the
   * one of the PointerMatrix.
   *
   * If <tt>M</tt> is zero, no
   * matrix is stored.
   */
  PointerMatrix (const MATRIX *M=0);

  /**
   * Constructor.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrix(const char *name);

  /**
   * Constructor. <tt>M</tt> points
   * to a matrix which must live
   * longer than the
   * PointerMatrix.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrix(const MATRIX *M,
                const char *name);

  // Use doc from base class
  virtual void clear();

  /**
   * Return whether the object is
   * empty.
   */
  bool empty () const;

  /**
   * Assign a new matrix
   * pointer. Deletes the old pointer
   * and releases its matrix.
   * @see SmartPointer
   */
  const PointerMatrix &operator= (const MATRIX *M);

  /**
   * Matrix-vector product.
   */
  virtual void vmult (VECTOR &dst,
                      const VECTOR &src) const;

  /**
   * Transposed matrix-vector product.
   */
  virtual void Tvmult (VECTOR &dst,
                       const VECTOR &src) const;

  /**
   * Matrix-vector product, adding to
   * <tt>dst</tt>.
   */
  virtual void vmult_add (VECTOR &dst,
                          const VECTOR &src) const;

  /**
   * Transposed matrix-vector product,
   * adding to <tt>dst</tt>.
   */
  virtual void Tvmult_add (VECTOR &dst,
                           const VECTOR &src) const;

private:
  /**
   * The pointer to the actual matrix.
   */
  SmartPointer<const MATRIX,PointerMatrix<MATRIX,VECTOR> > m;
};


/**
 * A pointer to be used as a matrix.  This class stores a pointer to a
 * matrix and can be used as a matrix itself in iterative methods.
 *
 * The main purpose for the existence of this class is its base class,
 * which only has a vector as template argument. Therefore, this
 * interface provides an abstract base class for matrices.
 *
 * This class differs form PointerMatrix by its additional
 * VectorMemory object and by the fact that it implements the
 * functions vmult_add() and Tvmult_add() only using vmult() and
 * Tvmult() of the MATRIX.
 *
 * @author Guido Kanschat 2006
 */
template<class MATRIX, class VECTOR>
class PointerMatrixAux : public PointerMatrixBase<VECTOR>
{
public:
  /**
   * Constructor.  The pointer in the
   * argument is stored in this
   * class. As usual, the lifetime of
   * <tt>*M</tt> must be longer than the
   * one of the PointerMatrixAux.
   *
   * If <tt>M</tt> is zero, no
   * matrix is stored.
   *
   * If <tt>mem</tt> is zero, then
   * GrowingVectorMemory
   * is used.
   */
  PointerMatrixAux (VectorMemory<VECTOR> *mem = 0,
                    const MATRIX *M=0);

  /**
   * Constructor not using a
   * matrix.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrixAux(VectorMemory<VECTOR> *mem,
                   const char *name);

  /**
   * Constructor. <tt>M</tt> points
   * to a matrix which must live
   * longer than the
   * PointerMatrixAux.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrixAux(VectorMemory<VECTOR> *mem,
                   const MATRIX *M,
                   const char *name);

  // Use doc from base class
  virtual void clear();

  /**
   * Return whether the object is
   * empty.
   */
  bool empty () const;

  /**
   * Assign a new VectorMemory
   * object for getting auxiliary
   * vectors.
   */
  void set_memory(VectorMemory<VECTOR> *mem);

  /**
   * Assign a new matrix
   * pointer. Deletes the old pointer
   * and releases its matrix.
   * @see SmartPointer
   */
  const PointerMatrixAux &operator= (const MATRIX *M);

  /**
   * Matrix-vector product.
   */
  virtual void vmult (VECTOR &dst,
                      const VECTOR &src) const;

  /**
   * Transposed matrix-vector product.
   */
  virtual void Tvmult (VECTOR &dst,
                       const VECTOR &src) const;

  /**
   * Matrix-vector product, adding to
   * <tt>dst</tt>.
   */
  virtual void vmult_add (VECTOR &dst,
                          const VECTOR &src) const;

  /**
   * Transposed matrix-vector product,
   * adding to <tt>dst</tt>.
   */
  virtual void Tvmult_add (VECTOR &dst,
                           const VECTOR &src) const;

private:
  /**
   * The backup memory if none was provided.
   */
  mutable GrowingVectorMemory<VECTOR> my_memory;

  /**
   * Object for getting the
   * auxiliary vector.
   */
  mutable SmartPointer<VectorMemory<VECTOR>,PointerMatrixAux<MATRIX,VECTOR> > mem;

  /**
   * The pointer to the actual matrix.
   */
  SmartPointer<const MATRIX,PointerMatrixAux<MATRIX,VECTOR> > m;
};



/**
 * Implement matrix multiplications for a vector using the
 * PointerMatrixBase functionality. Objects of this
 * class can be used in block matrices.
 *
 * Implements a matrix with image dimension 1 by using the scalar
 * product (#vmult()) and scalar multiplication (#Tvmult()) functions
 * of the Vector class.
 *
 * @author Guido Kanschat, 2006
 */
template <typename number>
class PointerMatrixVector : public PointerMatrixBase<Vector<number> >
{
public:
  /**
   * Constructor.  The pointer in the
   * argument is stored in this
   * class. As usual, the lifetime of
   * <tt>*M</tt> must be longer than the
   * one of the PointerMatrix.
   *
   * If <tt>M</tt> is zero, no
   * matrix is stored.
   */
  PointerMatrixVector (const Vector<number> *M=0);

  /**
   * Constructor.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrixVector (const char *name);

  /**
   * Constructor. <tt>M</tt> points
   * to a matrix which must live
   * longer than the
   * PointerMatrix.
   *
   * This class internally stores a pointer to a matrix via a SmartPointer
   * object. The SmartPointer class allows to associate a name with the
   * object pointed to that identifies the object that has the pointer,
   * in order to identify objects that still refer to the object pointed to.
   * The @p name argument to this function
   * is used to this end, i.e., you can in essence assign a name to
   * the current PointerMatrix object.
   */
  PointerMatrixVector (const Vector<number> *M,
                       const char *name);

  // Use doc from base class
  virtual void clear();

  /**
   * Return whether the object is
   * empty.
   */
  bool empty () const;

  /**
   * Assign a new matrix
   * pointer. Deletes the old pointer
   * and releases its matrix.
   * @see SmartPointer
   */
  const PointerMatrixVector &operator= (const Vector<number> *M);

  /**
   * Matrix-vector product,
   * actually the scalar product of
   * <tt>src</tt> and the vector
   * representing this matrix.
   *
   * The dimension of <tt>dst</tt>
   * is 1, while that of
   * <tt>src</tt> is the size of
   * the vector representing this
   * matrix.
   */
  virtual void vmult (Vector<number> &dst,
                      const Vector<number> &src) const;

  /**
   * Transposed matrix-vector
   * product, actually the
   * multiplication of the vector
   * representing this matrix with
   * <tt>src(0)</tt>.
   *
   * The dimension of <tt>drc</tt>
   * is 1, while that of
   * <tt>dst</tt> is the size of
   * the vector representing this
   * matrix.
   */
  virtual void Tvmult (Vector<number> &dst,
                       const Vector<number> &src) const;

  /**
   * Matrix-vector product, adding to
   * <tt>dst</tt>.
   *
   * The dimension of <tt>dst</tt>
   * is 1, while that of
   * <tt>src</tt> is the size of
   * the vector representing this
   * matrix.
   */
  virtual void vmult_add (Vector<number> &dst,
                          const Vector<number> &src) const;

  /**
   * Transposed matrix-vector product,
   * adding to <tt>dst</tt>.
   *
   * The dimension of <tt>src</tt>
   * is 1, while that of
   * <tt>dst</tt> is the size of
   * the vector representing this
   * matrix.
   */
  virtual void Tvmult_add (Vector<number> &dst,
                           const Vector<number> &src) const;

private:
  /**
   * The pointer to the actual matrix.
   */
  SmartPointer<const Vector<number>,PointerMatrixVector<number> > m;
};



/**
 * This function helps you creating a PointerMatrixBase object if you
 * do not want to provide the full template arguments of
 * PointerMatrix or PointerMatrixAux.
 *
 * Note that this function by default creates a PointerMatrixAux,
 * emulating the functions <tt>vmult_add</tt> and <tt>Tvmult_add</tt>,
 * using an auxiliary vector. It is overloaded for the library matrix
 * classes implementing these functions themselves. If you have such a
 * class, you should overload the function in order to save memory and
 * time.
 *
 * The result is a PointerMatrixBase* pointing to <tt>matrix</tt>. The
 * <TT>VECTOR</tt> argument is a dummy just used to determine the
 * template arguments.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrixAux
 */
template <class VECTOR, class MATRIX>
inline
PointerMatrixBase<VECTOR> *
new_pointer_matrix_base(MATRIX &matrix, const VECTOR &, const char *name = "PointerMatrixAux")
{
  return new PointerMatrixAux<MATRIX, VECTOR>(0, &matrix, name);
}

/**
 * Specialized version for IdentityMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const IdentityMatrix &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<IdentityMatrix, Vector<numberv> >(&matrix, name);
}


/**
 * Specialized version for FullMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const FullMatrix<numberm> &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<FullMatrix<numberm>, Vector<numberv> >(&matrix, name);
}


/**
 * Specialized version for LAPACKFullMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const LAPACKFullMatrix<numberm> &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<LAPACKFullMatrix<numberm>, Vector<numberv> >(&matrix, name);
}


/**
 * Specialized version for SparseMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const SparseMatrix<numberm> &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<SparseMatrix<numberm>, Vector<numberv> >(&matrix, name);
}


/**
 * Specialized version for BlockSparseMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <class VECTOR, typename numberm>
PointerMatrixBase<VECTOR> *
new_pointer_matrix_base(const BlockSparseMatrix<numberm> &matrix, const VECTOR &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<BlockSparseMatrix<numberm>, VECTOR>(&matrix, name);
}


/**
 * Specialized version for SparseMatrixEZ.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const SparseMatrixEZ<numberm> &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<SparseMatrixEZ<numberm>, Vector<numberv> >(&matrix, name);
}


/**
 * Specialized version for BlockSparseMatrixEZ.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <class VECTOR, typename numberm>
PointerMatrixBase<VECTOR> *
new_pointer_matrix_base(const BlockSparseMatrixEZ<numberm> &matrix, const VECTOR &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<BlockSparseMatrixEZ<numberm>, VECTOR>(&matrix, name);
}


/**
 * Specialized version for BlockMatrixArray.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<BlockVector<numberv> > *
new_pointer_matrix_base(const BlockMatrixArray<numberm> &matrix, const BlockVector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<BlockMatrixArray<numberm>, BlockVector<numberv> >(&matrix, name);
}


/**
 * Specialized version for TridiagonalMatrix.
 *
 * @relates PointerMatrixBase
 * @relates PointerMatrix
 */
template <typename numberv, typename numberm>
PointerMatrixBase<Vector<numberv> > *
new_pointer_matrix_base(const TridiagonalMatrix<numberm> &matrix, const Vector<numberv> &, const char *name = "PointerMatrix")
{
  return new PointerMatrix<TridiagonalMatrix<numberm>, Vector<numberv> >(&matrix, name);
}



/*@}*/
//---------------------------------------------------------------------------

template<class VECTOR>
inline
PointerMatrixBase<VECTOR>::~PointerMatrixBase ()
{}



//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (const MATRIX *M)
  : m(M, typeid(*this).name())
{}


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (const char *name)
  : m(0, name)
{}


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (
  const MATRIX *M,
  const char *name)
  : m(M, name)
{}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::clear ()
{
  m = 0;
}


template<class MATRIX, class VECTOR>
inline const PointerMatrix<MATRIX, VECTOR> &
PointerMatrix<MATRIX, VECTOR>::operator= (const MATRIX *M)
{
  m = M;
  return *this;
}


template<class MATRIX, class VECTOR>
inline bool
PointerMatrix<MATRIX, VECTOR>::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::vmult (VECTOR &dst,
                                      const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::Tvmult (VECTOR &dst,
                                       const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::vmult_add (VECTOR &dst,
                                          const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult_add (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::Tvmult_add (VECTOR &dst,
                                           const VECTOR &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult_add (dst, src);
}



//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR> *mem,
  const MATRIX *M)
  : mem(mem, typeid(*this).name()),
    m(M, typeid(*this).name())
{
  if (mem == 0) mem = &my_memory;
}


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR> *mem,
  const char *name)
  : mem(mem, name),
    m(0, name)
{
  if (mem == 0) mem = &my_memory;
}


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR> *mem,
  const MATRIX *M,
  const char *name)
  : mem(mem, name),
    m(M, name)
{
  if (mem == 0) mem = &my_memory;
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::clear ()
{
  m = 0;
}


template<class MATRIX, class VECTOR>
inline const PointerMatrixAux<MATRIX, VECTOR> &
PointerMatrixAux<MATRIX, VECTOR>::operator= (const MATRIX *M)
{
  m = M;
  return *this;
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::set_memory(VectorMemory<VECTOR> *M)
{
  mem = M;
  if (mem == 0)
    mem = &my_memory;
}


template<class MATRIX, class VECTOR>
inline bool
PointerMatrixAux<MATRIX, VECTOR>::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::vmult (VECTOR &dst,
                                         const VECTOR &src) const
{
  if (mem == 0)
    mem = &my_memory;
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::Tvmult (VECTOR &dst,
                                          const VECTOR &src) const
{
  if (mem == 0)
    mem = &my_memory;
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::vmult_add (VECTOR &dst,
                                             const VECTOR &src) const
{
  if (mem == 0)
    mem = &my_memory;
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  VECTOR *v = mem->alloc();
  v->reinit(dst);
  m->vmult (*v, src);
  dst += *v;
  mem->free(v);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::Tvmult_add (VECTOR &dst,
                                              const VECTOR &src) const
{
  if (mem == 0)
    mem = &my_memory;
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  VECTOR *v = mem->alloc();
  v->reinit(dst);
  m->Tvmult (*v, src);
  dst += *v;
  mem->free(v);
}


//----------------------------------------------------------------------//


template<typename number>
PointerMatrixVector<number>::PointerMatrixVector (const Vector<number> *M)
  : m(M, typeid(*this).name())
{}


template<typename number>
PointerMatrixVector<number>::PointerMatrixVector (const char *name)
  : m(0, name)
{}


template<typename number>
PointerMatrixVector<number>::PointerMatrixVector (
  const Vector<number> *M,
  const char *name)
  : m(M, name)
{}


template<typename number>
inline void
PointerMatrixVector<number>::clear ()
{
  m = 0;
}


template<typename number>
inline const PointerMatrixVector<number> &
PointerMatrixVector<number>::operator= (const Vector<number> *M)
{
  m = M;
  return *this;
}


template<typename number>
inline bool
PointerMatrixVector<number>::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<typename number>
inline void
PointerMatrixVector<number>::vmult (
  Vector<number> &dst,
  const Vector<number> &src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert (dst.size() == 1, ExcDimensionMismatch(dst.size(), 1));

  dst(0) = *m * src;
}


template<typename number>
inline void
PointerMatrixVector<number>::Tvmult (
  Vector<number> &dst,
  const Vector<number> &src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert(src.size() == 1, ExcDimensionMismatch(src.size(), 1));

  dst.equ (src(0), *m);
}


template<typename number>
inline void
PointerMatrixVector<number>::vmult_add (
  Vector<number> &dst,
  const Vector<number> &src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert (dst.size() == 1, ExcDimensionMismatch(dst.size(), 1));

  dst(0) += *m * src;
}


template<typename number>
inline void
PointerMatrixVector<number>::Tvmult_add (
  Vector<number> &dst,
  const Vector<number> &src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert(src.size() == 1, ExcDimensionMismatch(src.size(), 1));

  dst.add (src(0), *m);
}


DEAL_II_NAMESPACE_CLOSE

#endif
