//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
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

DEAL_II_NAMESPACE_OPEN

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
				      * Standard constructor. Matrices
				      * and the memory pool must be
				      * added later using
				      * initialize().
				      */
    ProductMatrix();
    
				     /**
				      * Constructor only assigning the
				      * memory pool. Matrices must be
				      * added by reinit() later.
				      */
    ProductMatrix(VectorMemory<VECTOR>& mem);
    
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
				      * Change the matrices and memory pool.
				      */
    template <class MATRIX1, class MATRIX2>
    void initialize(const MATRIX1& m1, const MATRIX2& m2,
		    VectorMemory<VECTOR>& mem);
    
				     /**
				      * Destructor.
				      */
    ~ProductMatrix();

				     // Doc in PointerMatrixBase
    void clear();
    
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
    SmartPointer<VectorMemory<VECTOR>,ProductMatrix<VECTOR> > mem;
				     /**
				      * Return some kind of
				      * identifier.
				      */
    virtual const void* get() const;
};


/**
 * A matrix that is the scaled version of another matrix.
 *
 * Matrix-vector products of this matrix are composed of those of the
 * original matrix and scaling by a constant factor.
 *
 * @author Guido Kanschat, 2007
 */
template<class VECTOR>
class ScaledMatrix : public Subscriptor
{
  public:
				     /**
				      * Constructor leaving an
				      * uninitialized object.
				      */
    ScaledMatrix ();
				     /**
				      * Constructor with initialization.
				      */
    template <class MATRIX>
    ScaledMatrix (const MATRIX& M, const double factor);

				     /**
				      * Destructor
				      */
    ~ScaledMatrix ();
				     /**
				      * Initialize for use with a new
				      * matrix and factor.
				      */
    template <class MATRIX>
    void initialize (const MATRIX& M, const double factor);

   				     /**
				      * Delete internal matrix pointer.
				      */
    void clear ();
    
				     /**
				      * Matrix-vector product.
				      */
    void vmult (VECTOR& w, const VECTOR& v) const;
    
				     /**
				      * Tranposed matrix-vector
				      * product.
				      */
    void Tvmult (VECTOR& w, const VECTOR& v) const;
    
  private:
				     /**
				      * The matrix.
				      */
    PointerMatrixBase<VECTOR>* m;
				     /**
				      * The scaling factor;
				      */
    double factor;
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
				      * Constructor leaving an
				      * unitialized
				      * matrix. initialize() must be
				      * called, before the matrix can
				      * be used.
				      */
    ProductSparseMatrix();
    
    void initialize(const MatrixType& m1,
		    const MatrixType& m2,
		    VectorMemory<VectorType>& mem);
    
				     // Doc in PointerMatrixBase
    void clear();
    
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
    SmartPointer<const MatrixType,ProductSparseMatrix<number,vector_number> > m1;
    
				     /**
				      * The right matrix of the product.
				      */
    SmartPointer<const MatrixType,ProductSparseMatrix<number,vector_number> > m2;
    
				     /**
				      * Memory for auxiliary vector.
				      */
    SmartPointer<VectorMemory<VectorType>,ProductSparseMatrix<number,vector_number>  > mem;
				     /**
				      * Return some kind of
				      * identifier.
				      */
    virtual const void* get() const;
};


/**
 * Mean value filter.  The vmult() functions of this matrix filter
 * out mean values of the vector.  If the vector is of type
 * BlockVector, then an additional parameter selects a single
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
    MeanValueFilter(unsigned int component = numbers::invalid_unsigned_int);

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
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 * 
 * @author Guido Kanschat, 2005
 */
template<class VECTOR>
class InverseMatrixRichardson : public Subscriptor
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
				      * Since we use two pointers, we
				      * must implement a destructor.
				      */
    ~InverseMatrixRichardson();
    
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
				      * A reference to the provided
				      * VectorMemory object.
				      */
    VectorMemory<VECTOR> &mem;
    
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
inline
ScaledMatrix<VECTOR>::ScaledMatrix()
		:
		m(0)
{}



template<class VECTOR>
template<class MATRIX>
inline
ScaledMatrix<VECTOR>::ScaledMatrix(const MATRIX& mat, const double factor)
		:
		m(new_pointer_matrix_base(mat, VECTOR())),
		factor(factor)
{}



template<class VECTOR>
template<class MATRIX>
inline
void
ScaledMatrix<VECTOR>::initialize(const MATRIX& mat, const double f)
{
  if (m) delete m;
  m = new_pointer_matrix_base(mat, VECTOR());
  factor = f;
}



template<class VECTOR>
inline
void
ScaledMatrix<VECTOR>::clear()
{
  if (m) delete m;
  m = 0;
}



template<class VECTOR>
inline
ScaledMatrix<VECTOR>::~ScaledMatrix()
{
  clear ();
}


template<class VECTOR>
inline
void
ScaledMatrix<VECTOR>::vmult (VECTOR& w, const VECTOR& v) const
{
  m->vmult(w, v);
  w *= factor;
}


template<class VECTOR>
inline
void
ScaledMatrix<VECTOR>::Tvmult (VECTOR& w, const VECTOR& v) const
{
  m->Tvmult(w, v);
  w *= factor;
}


//---------------------------------------------------------------------------

template<class VECTOR>
ProductMatrix<VECTOR>::ProductMatrix ()
  : m1(0), m2(0), mem(0)
{}


template<class VECTOR>
ProductMatrix<VECTOR>::ProductMatrix (VectorMemory<VECTOR>& m)
  : m1(0), m2(0), mem(&m)
{}


template<class VECTOR>
template<class MATRIX1, class MATRIX2>
ProductMatrix<VECTOR>::ProductMatrix (
  const MATRIX1& mat1,
  const MATRIX2& mat2,
  VectorMemory<VECTOR>& m)
  : mem(&m)
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
  if (m1) delete m1;
  if (m2) delete m2;
  m1 = new PointerMatrix<MATRIX1, VECTOR>(&mat1, typeid(*this).name());
  m2 = new PointerMatrix<MATRIX2, VECTOR>(&mat2, typeid(*this).name());
}


template<class VECTOR>
template<class MATRIX1, class MATRIX2>
void
ProductMatrix<VECTOR>::initialize (
  const MATRIX1& mat1,
  const MATRIX2& mat2,
  VectorMemory<VECTOR>& memory)
{
  mem = &memory;
  if (m1) delete m1;
  if (m2) delete m2;
  m1 = new PointerMatrix<MATRIX1, VECTOR>(&mat1, typeid(*this).name());
  m2 = new PointerMatrix<MATRIX2, VECTOR>(&mat2, typeid(*this).name());
}


template<class VECTOR>
ProductMatrix<VECTOR>::~ProductMatrix ()
{
  if (m1) delete m1;
  if (m2) delete m2;
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::clear ()
{
  if (m1) delete m1;
  m1 = 0;
  if (m2) delete m2;
  m2 = 0;
}


template<class VECTOR>
void
ProductMatrix<VECTOR>::vmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (mem != 0, ExcNotInitialized());
  Assert (m1 != 0, ExcNotInitialized());
  Assert (m2 != 0, ExcNotInitialized());
  
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
  Assert (mem != 0, ExcNotInitialized());
  Assert (m1 != 0, ExcNotInitialized());
  Assert (m2 != 0, ExcNotInitialized());
  
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
  Assert (mem != 0, ExcNotInitialized());
  Assert (m1 != 0, ExcNotInitialized());
  Assert (m2 != 0, ExcNotInitialized());
  
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
  Assert (mem != 0, ExcNotInitialized());
  Assert (m1 != 0, ExcNotInitialized());
  Assert (m2 != 0, ExcNotInitialized());
  
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
  if (matrix != 0)
    delete matrix;
  matrix = new PointerMatrix<MATRIX, VECTOR>(&m);
  if (precondition != 0)
    delete precondition;
  precondition = new PointerMatrix<PRECONDITION, VECTOR>(&p);;
}


DEAL_II_NAMESPACE_CLOSE

#endif
