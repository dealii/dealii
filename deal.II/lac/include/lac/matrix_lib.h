//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__matrix_lib_h
#define __deal2__matrix_lib_h

#include <base/subscriptor.h>
#include <lac/vector_memory.h>
#include <lac/pointer_matrix.h>
#include <lac/solver_richardson.h>

template<typename number> class Vector;
template<typename number> class BlockVector;
template<typename number> class SparseMatrix;

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
 * Here an example multiplying two different FullMatrix objects:
 * @include product_matrix.cc
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
				      * Change the matrices.
				      */
    template <class MATRIX1, class MATRIX2>
    void reinit(const MATRIX1& m1, const MATRIX2& m2);
    
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
				     /**
				      * Return some kind of
				      * identifier.
				      */
    virtual const void* get() const;
};


/**
 * Poor man's matrix product of two sparse matrices. Stores two
 * matrices #m1 and #m2 of arbitrary type SparseMatrix and implements
 * matrix-vector multiplications for the product
 * <i>M<sub>1</sub>M<sub>2</sub></i> by performing multiplication with
 * both factors consecutively.
 *
 * The documentation of ProductMatrix applies with exception that
 * these matrices here may be rectangular.
 *
 * @author Guido Kanschat, 2000, 2001, 2002, 2005
 */
template<typename number, typename vector_number>
class ProductSparseMatrix : public PointerMatrixBase<Vector<vector_number> >
{
  public:
				     /**
				      * Define the type of matrices used.
				      */
    typedef SparseMatrix<number> MatrixType;
    
				     /**
				      * Define the type of vectors we
				      * plly this matrix to.
				      */
    typedef Vector<vector_number> VectorType;
    
				     /**
				      * Constructor.  Additionally to
				      * the two constituting matrices, a
				      * memory pool for the auxiliary
				      * vector must be provided.
				      */
    ProductSparseMatrix(const MatrixType& m1,
			const MatrixType& m2,
			VectorMemory<VectorType>& mem);

				     /**
				      * Matrix-vector product <i>w =
				      * m1 * m2 * v</i>.
				      */
    virtual void vmult (VectorType&       w,
			const VectorType& v) const;
    
				     /**
				      * Tranposed matrix-vector
				      * product <i>w = m2<sup>T</sup>
				      * * m1<sup>T</sup> * v</i>.
				      */
    virtual void Tvmult (VectorType&       w,
			 const VectorType& v) const;
    
				     /**
				      * Adding matrix-vector product
				      * <i>w += m1 * m2 * v</i>
				      */
    virtual void vmult_add (VectorType&       w,
			    const VectorType& v) const;
    
				     /**
				      * Adding, tranposed
				      * matrix-vector product <i>w +=
				      * m2<sup>T</sup> *
				      * m1<sup>T</sup> * v</i>.
				      */
    virtual void Tvmult_add (VectorType&       w,
			     const VectorType& v) const;
    
  private:
				     /**
				      * The left matrix of the product.
				      */
    SmartPointer<const MatrixType> m1;
    
				     /**
				      * The right matrix of the product.
				      */
    SmartPointer<const MatrixType> m2;
    
				     /**
				      * Memory for auxiliary vector.
				      */
    SmartPointer<VectorMemory<VectorType> > mem;
				     /**
				      * Return some kind of
				      * identifier.
				      */
    virtual const void* get() const;
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



/**
 * Inverse matrix computed approximately by using the SolverRichardson
 * iterative solver. In particular, the function
 * SolverRichardson::Tsolve() allows for the implementation of
 * transpose matrix vector products.
 *
 * The functions vmult() and Tvmult() appoximate the inverse
 * iteratively starting with the vector <tt>dst</tt>. Functions
 * vmult_add() and Tvmult_add() start the iteration with a zero
 * vector.
 *
 * @ref Instantiations: some (Vector<float>, Vector<double>, BlockVector<float>, BlockVector<double>)
 * @author Guido Kanschat, 2005
 */
template<class VECTOR>
class InverseMatrixRichardson
{
  public:
				     /**
				      * Constructor, initializing the
				      * solver with a control and
				      * memory object. The inverted
				      * matrix and the preconditioner
				      * are added in initialize().
				      */
    InverseMatrixRichardson (SolverControl& control,
			     VectorMemory<VECTOR>& mem);
    
				     /**
				      * Initialization
				      * function. Provide a solver
				      * object, a matrix, and another
				      * preconditioner for this.
				      */
    template <class MATRIX, class PRECONDITION>
    void initialize (const MATRIX&,
		     const PRECONDITION&);

				     /**
				      * Access to the SolverControl
				      * object used by the solver.
				      */
    SolverControl& control() const;
				     /**
				      * Execute solver.
				      */
    void vmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Execute solver.
				      */
    void vmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * Execute transpose solver.
				      */
    void Tvmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Execute transpose solver.
				      */
    void Tvmult_add (VECTOR&, const VECTOR&) const;

  private:
				     /**
				      * Access to the provided
				      * VectorMemory object.
				      */
    mutable VectorMemory<VECTOR>& mem;
    
				     /**
				      * The solver object.
				      */
    mutable SolverRichardson<VECTOR> solver;
    
				     /**
				      * The matrix in use.
				      */
    PointerMatrixBase<VECTOR>* matrix;
    
				     /**
				      * The preconditioner to use.
				      */
    PointerMatrixBase<VECTOR>* precondition;
};




/*@}*/
//---------------------------------------------------------------------------

template<class VECTOR>
template<class MATRIX1, class MATRIX2>
ProductMatrix<VECTOR>::ProductMatrix (
  const MATRIX1& mat1,
  const MATRIX2& mat2,
  VectorMemory<VECTOR>& m)
  : mem(&m, typeid(*this).name())
{
  m1 = new PointerMatrix<MATRIX1, VECTOR>(&mat1, typeid(*this).name());
  m2 = new PointerMatrix<MATRIX2, VECTOR>(&mat2, typeid(*this).name());
}


template<class VECTOR>
template<class MATRIX1, class MATRIX2>
void
ProductMatrix<VECTOR>::reinit (
  const MATRIX1& mat1,
  const MATRIX2& mat2)
{
  delete m1;
  delete m2;
  m1 = new PointerMatrix<MATRIX1, VECTOR>(&mat1, typeid(*this).name());
  m2 = new PointerMatrix<MATRIX2, VECTOR>(&mat2, typeid(*this).name());
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


template<class VECTOR>
const void*
ProductMatrix<VECTOR>::get () const
{
  return (void*) m1;
}


//---------------------------------------------------------------------------

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

//-----------------------------------------------------------------------//

template <class VECTOR>
template <class MATRIX, class PRECONDITION>
inline void
InverseMatrixRichardson<VECTOR>::initialize (const MATRIX& m, const PRECONDITION& p)
{
  matrix = &m;
  precondition = &p;
}


#endif
