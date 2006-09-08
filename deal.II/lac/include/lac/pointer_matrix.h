//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__pointer_matrix_h
#define __deal2__pointer_matrix_h

#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/vector.h>

template<class VECTOR> class VectorMemory;

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
				      * the destructor of the derived
				      * class is called.
				      */
    virtual ~PointerMatrixBase ();

				     /**
				      * Reset pointer and release the
				      * matrix pointed to.
				      */
    virtual void clear () = 0;
    
				     /**
				      * Find out if two matrices point
				      * to the same object.
				      */
    bool operator == (const PointerMatrixBase<VECTOR>&) const;
    
				     /**
				      * Find out if two matrices do
				      * not point to the same object.
				      */
    bool operator != (const PointerMatrixBase<VECTOR>&) const;
    
				     /**
				      * Find out if this pointer is
				      * less.
				      */
    bool operator < (const PointerMatrixBase<VECTOR>&) const;
    
				     /**
				      * Find out if this pointer is
				      * less or equal.
				      */
    bool operator <= (const PointerMatrixBase<VECTOR>&) const;
    
				     /**
				      * Find out if this pointer is
				      * greater.
				      */
    bool operator > (const PointerMatrixBase<VECTOR>&) const;
    
				     /**
				      * Find out if this pointer is
				      * greater or equal.
				      */
    bool operator >= (const PointerMatrixBase<VECTOR>&) const;
    
    
				     /**
				      * Matrix-vector product.
				      */
    virtual void vmult (VECTOR& dst,
			const VECTOR& src) const = 0;
    
				     /**
				      * Tranposed matrix-vector product.
				      */
    virtual void Tvmult (VECTOR& dst,
			 const VECTOR& src) const = 0;
    
				     /**
				      * Matrix-vector product, adding to
				      * <tt>dst</tt>.
				      */
    virtual void vmult_add (VECTOR& dst,
			    const VECTOR& src) const = 0;
    
				     /**
				      * Tranposed matrix-vector product,
				      * adding to <tt>dst</tt>.
				      */
    virtual void Tvmult_add (VECTOR& dst,
			     const VECTOR& src) const = 0;

  private:
				     /**
				      * Get the pointer for comparison.
				      */
    virtual const void* get() const = 0;
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
    PointerMatrix (const MATRIX* M=0);

				     /**
				      * Constructor. The name argument
				      * is used to identify the
				      * SmartPointer for this object.
				      */
    PointerMatrix(const char* name);
    
				     /**
				      * Constructor. <tt>M</tt> points
				      * to a matrix which must live
				      * longer than the
				      * PointerMatrix. The name
				      * argument is used to identify
				      * the SmartPointer for this
				      * object.
				      */
    PointerMatrix(const MATRIX* M,
		  const char* name);

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
    const PointerMatrix& operator= (const MATRIX* M);
    
				     /**
				      * Matrix-vector product.
				      */
    virtual void vmult (VECTOR& dst,
			const VECTOR& src) const;
    
				     /**
				      * Tranposed matrix-vector product.
				      */
    virtual void Tvmult (VECTOR& dst,
			 const VECTOR& src) const;
    
				     /**
				      * Matrix-vector product, adding to
				      * <tt>dst</tt>.
				      */
    virtual void vmult_add (VECTOR& dst,
			    const VECTOR& src) const;
    
				     /**
				      * Tranposed matrix-vector product,
				      * adding to <tt>dst</tt>.
				      */
    virtual void Tvmult_add (VECTOR& dst,
			     const VECTOR& src) const;
    
  private:
				     /**
				      * Return the address of the
				      * matrix for comparison.
				      */
    virtual const void* get() const;

				     /**
				      * The pointer to the actual matrix.
				      */
    SmartPointer<const MATRIX> m;
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
				      */
    PointerMatrixAux (VectorMemory<VECTOR>* mem = 0,
		      const MATRIX* M=0);

				     /**
				      * Constructor. The name argument
				      * is used to identify the
				      * SmartPointer for this object.
				      */
    PointerMatrixAux(VectorMemory<VECTOR>* mem,
		     const char* name);
    
				     /**
				      * Constructor. <tt>M</tt> points
				      * to a matrix which must live
				      * longer than the
				      * PointerMatrixAux. The name
				      * argument is used to identify
				      * the SmartPointer for this
				      * object.
				      */
    PointerMatrixAux(VectorMemory<VECTOR>* mem,
		     const MATRIX* M,
		     const char* name);

				     // Use doc from base class
    virtual void clear();
    
				     /**
				      * Return whether the object is
				      * empty. 
				      */
    bool empty () const;

				     /**
				      * Assign a new VectorMemory
				     object for getting auxiliary vectors.
				      */
    void set_memory(VectorMemory<VECTOR>* mem);
    
				     /**
				      * Assign a new matrix
				      * pointer. Deletes the old pointer
				      * and releases its matrix.
				      * @see SmartPointer
				      */
    const PointerMatrixAux& operator= (const MATRIX* M);
    
				     /**
				      * Matrix-vector product.
				      */
    virtual void vmult (VECTOR& dst,
			const VECTOR& src) const;
    
				     /**
				      * Tranposed matrix-vector product.
				      */
    virtual void Tvmult (VECTOR& dst,
			 const VECTOR& src) const;
    
				     /**
				      * Matrix-vector product, adding to
				      * <tt>dst</tt>.
				      */
    virtual void vmult_add (VECTOR& dst,
			    const VECTOR& src) const;
    
				     /**
				      * Tranposed matrix-vector product,
				      * adding to <tt>dst</tt>.
				      */
    virtual void Tvmult_add (VECTOR& dst,
			     const VECTOR& src) const;
    
  private:
				     /**
				      * Return the address of the
				      * matrix for comparison.
				      */
    virtual const void* get() const;

				     /**
				      * Object for getting the
				      * auxiliary vector.
				      */
    SmartPointer<VectorMemory<VECTOR> > mem;
    
				     /**
				      * The pointer to the actual matrix.
				      */
    SmartPointer<const MATRIX> m;
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
 * @author Guidl Kanschat, 2006
 */
template <typename number>
PointerMatrixVector : public PointerMatrixBase<Vector<number> >
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
    PointerMatrix (const Vector<number>* M=0);

				     /**
				      * Constructor. The name argument
				      * is used to identify the
				      * SmartPointer for this object.
				      */
    PointerMatrix(const char* name);
    
				     /**
				      * Constructor. <tt>M</tt> points
				      * to a matrix which must live
				      * longer than the
				      * PointerMatrix. The name
				      * argument is used to identify
				      * the SmartPointer for this
				      * object.
				      */
    PointerMatrix(const Vector<number>* M,
		  const char* name);

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
    const PointerMatrix& operator= (const Vector<number>* M);
    
				     /**
				      * Matrix-vector product,
				      * actually the scalar product of
				      * <tt>src</tt> and #m.
				      */
    virtual void vmult (Vector<number>& dst,
			const Vector<number>& src) const;
    
				     /**
				      * Tranposed matrix-vector
				      * product, actually the
				      * multiplication of #m with
				      * <tt>src(0)</tt>.
				      */
    virtual void Tvmult (Vector<number>& dst,
			 const Vector<number>& src) const;
    
				     /**
				      * Matrix-vector product, adding to
				      * <tt>dst</tt>.
				      */
    virtual void vmult_add (Vector<number>& dst,
			    const Vector<number>& src) const;
    
				     /**
				      * Tranposed matrix-vector product,
				      * adding to <tt>dst</tt>.
				      */
    virtual void Tvmult_add (Vector<number>& dst,
			     const Vector<number>& src) const;
    
  private:
				     /**
				      * Return the address of the
				      * matrix for comparison.
				      */
    virtual const void* get() const;

				     /**
				      * The pointer to the actual matrix.
				      */
    SmartPointer<const Vector<number> > m;
};



/**
 * This function helps you creating a PointerMatrixBase object if you
 * do not want to provide the full template arguments of
 * PointerMatrix.
 *
 * The result is a PointerMatrixBase* pointing to <tt>matrix</tt>. The
 * <TT>VECTOR</tt> argument is a dummy just used to determine the
 * template arguments.
 *
 * @relates PointerMatrixBase
 */
template <class MATRIX, class VECTOR>
inline
PointerMatrixBase<VECTOR>*
new_pointer_matrix_base(MATRIX& matrix, const VECTOR&)
{
  return new PointerMatrix<MATRIX, VECTOR>(&matrix);
}

/*@}*/
//---------------------------------------------------------------------------

template<class VECTOR>
inline 
PointerMatrixBase<VECTOR>::~PointerMatrixBase ()
{}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator == (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() == other.get());
}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator != (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() != other.get());
}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator < (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() < other.get());
}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator <= (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() <= other.get());
}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator > (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() > other.get());
}



template<class VECTOR>
inline
bool
PointerMatrixBase<VECTOR>::operator >= (const PointerMatrixBase<VECTOR>& other) const
{
  return (get() >= other.get());
}


//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (const MATRIX* M)
  : m(M, typeid(*this).name())
{}


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (const char* name)
  : m(0, name)
{}


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (
  const MATRIX* M,
  const char* name)
  : m(M, name)
{}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::clear ()
{
  m = 0;
}


template<class MATRIX, class VECTOR>
inline const PointerMatrix<MATRIX, VECTOR>&
PointerMatrix<MATRIX, VECTOR>::operator= (const MATRIX* M)
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
PointerMatrix<MATRIX, VECTOR>::vmult (VECTOR& dst,
				      const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::Tvmult (VECTOR& dst,
				       const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::vmult_add (VECTOR& dst,
					  const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult_add (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrix<MATRIX, VECTOR>::Tvmult_add (VECTOR& dst,
					   const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult_add (dst, src);
}


template<class MATRIX, class VECTOR>
inline const void*
PointerMatrix<MATRIX, VECTOR>::get () const
{
  return m;
}


//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR>* mem,
  const MATRIX* M)
		: mem(mem, typeid(*this).name()),
		  m(M, typeid(*this).name())
{}


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR>* mem,
  const char* name)
		: mem(mem, name),
		  m(0, name)
{}


template<class MATRIX, class VECTOR>
PointerMatrixAux<MATRIX, VECTOR>::PointerMatrixAux (
  VectorMemory<VECTOR>* mem,
  const MATRIX* M,
  const char* name)
		: mem(mem, name),
		  m(M, name)
{}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::clear ()
{
  m = 0;
}


template<class MATRIX, class VECTOR>
inline const PointerMatrixAux<MATRIX, VECTOR>&
PointerMatrixAux<MATRIX, VECTOR>::operator= (const MATRIX* M)
{
  m = M;
  return *this;
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::set_memory(VectorMemory<VECTOR>* M)
{
  mem = M;
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
PointerMatrixAux<MATRIX, VECTOR>::vmult (VECTOR& dst,
				      const VECTOR& src) const
{
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::Tvmult (VECTOR& dst,
				       const VECTOR& src) const
{
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::vmult_add (VECTOR& dst,
					  const VECTOR& src) const
{
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m->vmult (*v, src);
  dst += *v;
  mem->free(v);
}


template<class MATRIX, class VECTOR>
inline void
PointerMatrixAux<MATRIX, VECTOR>::Tvmult_add (VECTOR& dst,
					      const VECTOR& src) const
{
  Assert (mem != 0, ExcNotInitialized());
  Assert (m != 0, ExcNotInitialized());
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m->Tvmult (*v, src);
  dst += *v;
  mem->free(v);
}


template<class MATRIX, class VECTOR>
inline const void*
PointerMatrixAux<MATRIX, VECTOR>::get () const
{
  return m;
}


//----------------------------------------------------------------------//


template<typename number>
PointerMatrixVector<Vector<number> >::PointerMatrix (const MATRIX* M)
  : m(M, typeid(*this).name())
{}


template<typename number>
PointerMatrixVector<Vector<number> >::PointerMatrix (const char* name)
  : m(0, name)
{}


template<typename number>
PointerMatrixVector<Vector<number> >::PointerMatrix (
  const MATRIX* M,
  const char* name)
  : m(M, name)
{}


template<typename number>
inline void
PointerMatrixVector<Vector<number> >::clear ()
{
  m = 0;
}


template<typename number>
inline const PointerMatrixVector<Vector<number> >&
PointerMatrixVector<Vector<number> >::operator= (const MATRIX* M)
{
  m = M;
  return *this;
}


template<typename number>
inline bool
PointerMatrixVector<Vector<number> >::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<typename number>
inline void
PointerMatrixVector<Vector<number> >::vmult (
  Vector<number>& dst,
  const Vector<number>& src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert (dst.size() == 1, ExcDimensionMismatch(dst.size(), 1));
  
  dst(0) = *m * src;
}


template<typename number>
inline void
PointerMatrixVector<Vector<number> >::Tvmult (
  Vector<number>& dst,
  const Vector<number>& src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert(src.size() == 1, ExcDimensionMismatch(src.size(), 1));
  
  dst.equ (src(0), *m);
}


template<typename number>
inline void
PointerMatrixVector<Vector<number> >::vmult_add (
  Vector<number>& dst,
  const Vector<number>& src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert (dst.size() == 1, ExcDimensionMismatch(dst.size(), 1));
  
  dst(0) += *m * src;
}


template<typename number>
inline void
PointerMatrixVector<Vector<number> >::Tvmult_add (
  Vector<number>& dst,
  const Vector<number>& src) const
{
  Assert (m != 0, ExcNotInitialized());
  Assert(src.size() == 1, ExcDimensionMismatch(src.size(), 1));
  
  dst.add (src(0), *m);
}


template<typename number>
inline const void*
PointerMatrixVector<Vector<number> >::get () const
{
  return m;
}



#endif
