// $Id$

#ifndef __deal2__matrix_lib_h
#define __deal2__matrix_lib_h

#include <base/subscriptor.h>
#include <lac/vector_memory.h>
#include <lac/pointer_matrix.h>

template<typename number> class Vector;
template<typename number> class BlockVector;

/*! @addtogroup Matrix2
 *@{
 */


/**
 * Poor man's matrix product of two quadratic matrices. Stores two
 * quadratic matrices #m1 and #m2 of arbitrary types and implements
 * matrix-vector multiplications for the product
 * <i>M<sub>1</sub>M<sub>2</sub></i> by performing multiplication with
 * both factors consecutively.
 *
 * @author Guido Kanschat, 2000, 2001, 2002, 2005
 */
template<class VECTOR>
class ProductMatrix : public PointerMatrixBase<VECTOR>
{
  public:
				     /**
				      * Constructor.  Additionally to
				      * the two constituting matrices, a
				      * memory pool for the auxiliary
				      * vector must be provided.
				      */
    template <class MATRIX1, class MATRIX2>
    ProductMatrix(const MATRIX1& m1,
		  const MATRIX2& m2,
		  VectorMemory<VECTOR>& mem);

				     /**
				      * Destructor.
				      */
    ~ProductMatrix();
    
				     /**
				      * Matrix-vector product <i>w =
				      * m1 * m2 * v</i>.
				      */
    virtual void vmult (VECTOR&       w,
			const VECTOR& v) const;
    
				     /**
				      * Tranposed matrix-vector
				      * product <i>w = m2<sup>T</sup>
				      * * m1<sup>T</sup> * v</i>.
				      */
    virtual void Tvmult (VECTOR&       w,
			 const VECTOR& v) const;
    
				     /**
				      * Adding matrix-vector product
				      * <i>w += m1 * m2 * v</i>
				      */
    virtual void vmult_add (VECTOR&       w,
			    const VECTOR& v) const;
    
				     /**
				      * Adding, tranposed
				      * matrix-vector product <i>w +=
				      * m2<sup>T</sup> *
				      * m1<sup>T</sup> * v</i>.
				      */
    virtual void Tvmult_add (VECTOR&       w,
			     const VECTOR& v) const;
    
  private:
				     /**
				      * The left matrix of the product.
				      */
    PointerMatrixBase<VECTOR>* m1;
    
				     /**
				      * The right matrix of the product.
				      */
    PointerMatrixBase<VECTOR>* m2;
    
				     /**
				      * Memory for auxiliary vector.
				      */
    SmartPointer<VectorMemory<VECTOR> > mem;
};


/**
 * Mean value filter.  The @p vmult functions of this matrix filter
 * out mean values of the vector.  If the vector is of type
 * @p BlockVector, then an additional parameter selects a single
 * component for this operation.
 *
 * @author Guido Kanschat, 2002, 2003
 */
class MeanValueFilter : public Subscriptor
{
  public:
				     /**
				      * Constructor, optionally
				      * selecting a component.
				      */
    MeanValueFilter(unsigned int component = deal_II_numbers::invalid_unsigned_int);

				     /**
				      * Subtract mean value from @p v.
				      */
    template <typename number>
    void filter (Vector<number>& v) const;
    
				     /**
				      * Subtract mean value from @p v.
				      */
    template <typename number>
    void filter (BlockVector<number>& v) const;
    
				     /**
				      * Return the source vector with
				      * subtracted mean value.
				      */
    template <typename number>
    void vmult (Vector<number>& dst,
		const Vector<number>& src) const;
    
				     /**
				      * Add source vector with
				      * subtracted mean value to dest.
				      */
    template <typename number>
    void vmult_add (Vector<number>& dst,
		    const Vector<number>& src) const;
    
				     /**
				      * Return the source vector with
				      * subtracted mean value in
				      * selected component.
				      */
    template <typename number>
    void vmult (BlockVector<number>& dst,
		const BlockVector<number>& src) const;
    
				     /**
				      * Add a soruce to dest, where
				      * the mean value in the selected
				      * component is subtracted.
				      */
    template <typename number>
    void vmult_add (BlockVector<number>& dst,
		    const BlockVector<number>& src) const;


				     /**
				      * Not implemented.
				      */
    template <typename VECTOR>
    void Tvmult(VECTOR&, const VECTOR&) const;
    
				     /**
				      * Not implemented.
				      */
    template <typename VECTOR>
    void Tvmult_add(VECTOR&, const VECTOR&) const;
    
  private:
				     /**
				      * Component for filtering block vectors.
				      */
    unsigned int component;
};

/*@}*/
//----------------------------------------------------------------------//

template<class VECTOR>
template<class MATRIX1, class MATRIX2>
ProductMatrix<VECTOR>::ProductMatrix (
  const MATRIX1& mat1,
  const MATRIX2& mat2,
  VectorMemory<VECTOR>& m)
  : mem(&m)
{
  m1 = new PointerMatrix<MATRIX1, VECTOR>(&mat1);
  m2 = new PointerMatrix<MATRIX2, VECTOR>(&mat2);
}


template<class VECTOR>
ProductMatrix<VECTOR>::~ProductMatrix ()
{
  delete m1;
  delete m2;
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::vmult (VECTOR& dst, const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m2->vmult (*v, src);
  m1->vmult (dst, *v);
  mem->free(v);
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::vmult_add (VECTOR& dst, const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m2->vmult (*v, src);
  m1->vmult_add (dst, *v);
  mem->free(v);
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::Tvmult (VECTOR& dst, const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m1->Tvmult (*v, src);
  m2->Tvmult (dst, *v);
  mem->free(v);
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::Tvmult_add (VECTOR& dst, const VECTOR& src) const
{
  VECTOR* v = mem->alloc();
  v->reinit(dst);
  m1->Tvmult (*v, src);
  m2->Tvmult_add (dst, *v);
  mem->free(v);
}


//----------------------------------------------------------------------//

template <class VECTOR>
inline void
MeanValueFilter::Tvmult(VECTOR&, const VECTOR&) const
{
  Assert(false, ExcNotImplemented());
}


template <class VECTOR>
inline void
MeanValueFilter::Tvmult_add(VECTOR&, const VECTOR&) const
{
  Assert(false, ExcNotImplemented());
}


#endif
