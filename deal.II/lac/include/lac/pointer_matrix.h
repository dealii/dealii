// $Id$

#ifndef __deal2__pointer_matrix_h
#define __deal2__pointer_matrix_h

#include <base/subscriptor.h>
#include <lac/vector_memory.h>

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
				    * Virtual destructor.  Does
				    * nothing except making sure that
				    * the destructor of the derived
				    * class is called.
				    */
  virtual ~PointerMatrixBase ();

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
				    * @p{dst}.
				    */
  virtual void vmult_add (VECTOR& dst,
			  const VECTOR& src) const = 0;

				   /**
				    * Tranposed matrix-vector product,
				    * adding to @p{dst}.
				    */
  virtual void Tvmult_add (VECTOR& dst,
			   const VECTOR& src) const = 0;
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
				    * @p{*M} must be longer than the
				    * one of the @p{PointerMatrix}.
				    *
				    * If @p{M} is zero, no matrix is stored.
				    */
  PointerMatrix (const MATRIX* M=0);

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
				    * @p{dst}.
				    */
  virtual void vmult_add (VECTOR& dst,
			  const VECTOR& src) const;
  
				   /**
				    * Tranposed matrix-vector product,
				    * adding to @p{dst}.
				    */
  virtual void Tvmult_add (VECTOR& dst,
			   const VECTOR& src) const;
  
private:
				   /**
				    * The pointer to the actual matrix.
				    */
  SmartPointer<const MATRIX> m;
};


/**
 * Poor man's matrix product. Stores two matrices $m_1$ and $m_2$ and
 * implements matrix-vector multiplications for the product $m_1 m_2$
 * by performing multiplication with both factors consecutively.
 *
 * @author Guido Kanschat, 2000, 2001, 2002
 */
template<class MATRIX1, class MATRIX2, class VECTOR>
class ProductMatrix : public PointerMatrixBase<VECTOR>
{
  public:
				   /**
				    * Constructor.  Additionally to
				    * the two constituting matrices, a
				    * memory pool for the auxiliary
				    * vector must be provided.
				    */
  ProductMatrix(const MATRIX1& m1,
		const MATRIX2& m2,
		VectorMemory<VECTOR>& mem);
  
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
				    * @p{dst}.
				    */
  virtual void vmult_add (VECTOR& dst,
			  const VECTOR& src) const;
  
				   /**
				    * Tranposed matrix-vector product,
				    * adding to @p{dst}.
				    */
  virtual void Tvmult_add (VECTOR& dst,
			   const VECTOR& src) const;
  
  private:
				   /**
				    * The left matrix of the product.
				    */
  SmartPointer<const MATRIX1> m1;

				   /**
				    * The right matrix of the product.
				    */
  SmartPointer<const MATRIX2> m2;

				   /**
				    * Memory for auxiliary vector.
				    */
  SmartPointer<VectorMemory<VECTOR> > mem;
};


//----------------------------------------------------------------------//

template<class VECTOR>
inline 
PointerMatrixBase<VECTOR>::~PointerMatrixBase ()
{}


template<class MATRIX, class VECTOR>
PointerMatrix<MATRIX, VECTOR>::PointerMatrix (const MATRIX* M)
  :
  m(M)
{}

template<class MATRIX, class VECTOR>
const PointerMatrix<MATRIX, VECTOR>&
PointerMatrix<MATRIX, VECTOR>::operator= (const MATRIX* M)
{
  m = M;
  return *this;
}

template<class MATRIX, class VECTOR>
void
PointerMatrix<MATRIX, VECTOR>::vmult (VECTOR& dst,
				      const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<class MATRIX, class VECTOR>
void
PointerMatrix<MATRIX, VECTOR>::Tvmult (VECTOR& dst,
				       const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<class MATRIX, class VECTOR>
void
PointerMatrix<MATRIX, VECTOR>::vmult_add (VECTOR& dst,
					  const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult_add (dst, src);
}


template<class MATRIX, class VECTOR>
void
PointerMatrix<MATRIX, VECTOR>::Tvmult_add (VECTOR& dst,
					   const VECTOR& src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult_add (dst, src);
}


//----------------------------------------------------------------------//

template<class MATRIX1, class MATRIX2, class VECTOR>
ProductMatrix<MATRIX1, MATRIX2, VECTOR>::ProductMatrix (
  const MATRIX1& mat1,
  const MATRIX2& mat2,
  VectorMemory<VECTOR>& m)
  : m1(&mat1),
    m2(&mat2),
    mem(&m)
{}


template<class MATRIX1, class MATRIX2, class VECTOR>
void
ProductMatrix<MATRIX1, MATRIX2, VECTOR>::vmult (VECTOR& dst,
				      const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m2->vmult (*v, src);
  m1->vmult (dst, *v);
  mem->free(v);
}


template<class MATRIX1, class MATRIX2, class VECTOR>
void
ProductMatrix<MATRIX1, MATRIX2, VECTOR>::vmult_add (VECTOR& dst,
					  const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m2->vmult (*v, src);
  m1->vmult_add (dst, *v);
  mem->free(v);
}


template<class MATRIX1, class MATRIX2, class VECTOR>
void
ProductMatrix<MATRIX1, MATRIX2, VECTOR>::Tvmult (
  VECTOR& dst,
  const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m1->Tvmult (*v, src);
  m2->Tvmult (dst, *v);
  mem->free(v);
}


template<class MATRIX1, class MATRIX2, class VECTOR>
void
ProductMatrix<MATRIX1, MATRIX2, VECTOR>::Tvmult_add (
  VECTOR& dst,
  const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m1->Tvmult (*v, src);
  m2->Tvmult_add (dst, *v);
  mem->free(v);
}


#endif
