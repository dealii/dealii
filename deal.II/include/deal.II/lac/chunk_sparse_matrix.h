//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name:  $
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__chunk_sparse_matrix_h
#define __deal2__chunk_sparse_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/chunk_sparsity_pattern.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

template<typename number> class Vector;
template<typename number> class FullMatrix;

/*! @addtogroup Matrix1
 *@{
 */


/**
 * Sparse matrix. This class implements the function to store values
 * in the locations of a sparse matrix denoted by a
 * SparsityPattern. The separation of sparsity pattern and values is
 * done since one can store data elements of different type in these
 * locations without the SparsityPattern having to know this, and more
 * importantly one can associate more than one matrix with the same
 * sparsity pattern.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Wolfgang Bangerth, 2008
 */
template <typename number>
class ChunkSparseMatrix : public virtual Subscriptor
{
  public:
				     /**
				      * Type of matrix entries. In analogy to
				      * the STL container classes.
				      */
    typedef number value_type;

				     /**
				      * Declare a type that has holds
				      * real-valued numbers with the
				      * same precision as the template
				      * argument to this class. If the
				      * template argument of this
				      * class is a real data type,
				      * then real_type equals the
				      * template argument. If the
				      * template argument is a
				      * std::complex type then
				      * real_type equals the type
				      * underlying the complex
				      * numbers.
				      *
				      * This typedef is used to
				      * represent the return type of
				      * norms.
				      */
    typedef typename numbers::NumberTraits<number>::real_type real_type;    

                                     /**
                                      * A structure that describes some of the
                                      * traits of this class in terms of its
                                      * run-time behavior. Some other classes
                                      * (such as the block matrix classes)
                                      * that take one or other of the matrix
                                      * classes as its template parameters can
                                      * tune their behavior based on the
                                      * variables in this class.
                                      */
    struct Traits
    {
                                         /**
                                          * It is safe to elide additions of
                                          * zeros to individual elements of
                                          * this matrix.
                                          */
        static const bool zero_addition_can_be_elided = true;
    };
    
/**
 * @name Constructors and initalization.
 */
//@{
				     /**
				      * Constructor; initializes the matrix to
				      * be empty, without any structure, i.e.
				      * the matrix is not usable at all. This
				      * constructor is therefore only useful
				      * for matrices which are members of a
				      * class. All other matrices should be
				      * created at a point in the data flow
				      * where all necessary information is
				      * available.
				      *
				      * You have to initialize
				      * the matrix before usage with
				      * reinit(const ChunkSparsityPattern&).
				      */
    ChunkSparseMatrix ();

				     /**
				      * Copy constructor. This constructor is
				      * only allowed to be called if the matrix
				      * to be copied is empty. This is for the
				      * same reason as for the
				      * ChunkSparsityPattern, see there for the
				      * details.
				      *
				      * If you really want to copy a whole
				      * matrix, you can do so by using the
				      * copy_from() function.
				      */
    ChunkSparseMatrix (const ChunkSparseMatrix &);

				     /**
				      * Constructor. Takes the given
				      * matrix sparsity structure to
				      * represent the sparsity pattern
				      * of this matrix. You can change
				      * the sparsity pattern later on
				      * by calling the reinit(const
				      * ChunkSparsityPattern&) function.
				      *
				      * You have to make sure that the
				      * lifetime of the sparsity
				      * structure is at least as long
				      * as that of this matrix or as
				      * long as reinit(const
				      * ChunkSparsityPattern&) is not
				      * called with a new sparsity
				      * pattern.
				      *
				      * The constructor is marked
				      * explicit so as to disallow
				      * that someone passes a sparsity
				      * pattern in place of a sparse
				      * matrix to some function, where
				      * an empty matrix would be
				      * generated then.
				      */
    explicit ChunkSparseMatrix (const ChunkSparsityPattern &sparsity);

				     /**
				      * Copy constructor: initialize
				      * the matrix with the identity
				      * matrix. This constructor will
				      * throw an exception if the
				      * sizes of the sparsity pattern
				      * and the identity matrix do not
				      * coincide, or if the sparsity
				      * pattern does not provide for
				      * nonzero entries on the entire
				      * diagonal.
				      */
    ChunkSparseMatrix (const ChunkSparsityPattern &sparsity,
		       const IdentityMatrix  &id);
    
				     /**
				      * Destructor. Free all memory, but do not
				      * release the memory of the sparsity
				      * structure.
				      */
    virtual ~ChunkSparseMatrix ();

				     /** 
				      * Copy operator. Since copying
				      * entire sparse matrices is a
				      * very expensive operation, we
				      * disallow doing so except for
				      * the special case of empty
				      * matrices of size zero. This
				      * doesn't seem particularly
				      * useful, but is exactly what
				      * one needs if one wanted to
				      * have a
				      * <code>std::vector@<ChunkSparseMatrix@<double@>
				      * @></code>: in that case, one
				      * can create a vector (which
				      * needs the ability to copy
				      * objects) of empty matrices
				      * that are then later filled
				      * with something useful.
				      */
    ChunkSparseMatrix<number>& operator = (const ChunkSparseMatrix<number> &);

				     /**
				      * Copy operator: initialize
				      * the matrix with the identity
				      * matrix. This operator will
				      * throw an exception if the
				      * sizes of the sparsity pattern
				      * and the identity matrix do not
				      * coincide, or if the sparsity
				      * pattern does not provide for
				      * nonzero entries on the entire
				      * diagonal.
				      */
    ChunkSparseMatrix<number> &
    operator= (const IdentityMatrix  &id);
    
                                     /**
                                      * This operator assigns a scalar to
                                      * a matrix. Since this does usually
                                      * not make much sense (should we set
                                      * all matrix entries to this value?
                                      * Only the nonzero entries of the
                                      * sparsity pattern?), this operation
                                      * is only allowed if the actual
                                      * value to be assigned is zero. This
                                      * operator only exists to allow for
                                      * the obvious notation
                                      * <tt>matrix=0</tt>, which sets all
                                      * elements of the matrix to zero,
                                      * but keep the sparsity pattern
                                      * previously used.
                                      */
    ChunkSparseMatrix & operator = (const double d);

				     /**
				      * Reinitialize the sparse matrix
				      * with the given sparsity
				      * pattern. The latter tells the
				      * matrix how many nonzero
				      * elements there need to be
				      * reserved.
				      *
				      * Regarding memory allocation,
				      * the same applies as said
				      * above.
				      *
				      * You have to make sure that the
				      * lifetime of the sparsity
				      * structure is at least as long
				      * as that of this matrix or as
				      * long as reinit(const
				      * ChunkSparsityPattern &) is not
				      * called with a new sparsity
				      * structure.
				      *
				      * The elements of the matrix are
				      * set to zero by this function.
				      */
    virtual void reinit (const ChunkSparsityPattern &sparsity);

				     /**
				      * Release all memory and return
				      * to a state just like after
				      * having called the default
				      * constructor. It also forgets
				      * the sparsity pattern it was
				      * previously tied to.
				      */
    virtual void clear ();
//@}    
/**
 * @name Information on the matrix
 */
//@{
				     /**
				      * Return whether the object is
				      * empty. It is empty if either
				      * both dimensions are zero or no
				      * ChunkSparsityPattern is
				      * associated.
				      */
    bool empty () const;

				     /**
				      * Return the dimension of the
				      * image space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int m () const;
    
				     /**
				      * Return the dimension of the
				      * range space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int n () const;

				     /**
				      * Return the number of nonzero
				      * elements of this
				      * matrix. Actually, it returns
				      * the number of entries in the
				      * sparsity pattern; if any of
				      * the entries should happen to
				      * be zero, it is counted anyway.
				      */
    unsigned int n_nonzero_elements () const;

				     /**
				      * Return the number of actually
				      * nonzero elements of this
				      * matrix.
				      *
				      * Note, that this function does
				      * (in contrary to
				      * n_nonzero_elements()) not
				      * count all entries of the
				      * sparsity pattern but only the
				      * ones that are nonzero.
				      */
    unsigned int n_actually_nonzero_elements () const;
    
				     /**
				      * Return a (constant) reference
				      * to the underlying sparsity
				      * pattern of this matrix.
				      *
				      * Though the return value is
				      * declared <tt>const</tt>, you
				      * should be aware that it may
				      * change if you call any
				      * nonconstant function of
				      * objects which operate on it.
				      */
    const ChunkSparsityPattern & get_sparsity_pattern () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object. See
				      * MemoryConsumption.
				      */
    std::size_t memory_consumption () const;
    
//@}
/**
 * @name Modifying entries
 */
//@{
				     /**
				      * Set the element (<i>i,j</i>)
				      * to <tt>value</tt>. Throws an
				      * error if the entry does not
				      * exist or if <tt>value</tt> is
				      * not a finite number. Still, it
				      * is allowed to store zero
				      * values in non-existent fields.
				      */
    void set (const unsigned int i,
              const unsigned int j,
	      const number value);
    
				     /**
				      * Add <tt>value</tt> to the
				      * element (<i>i,j</i>).  Throws
				      * an error if the entry does not
				      * exist or if <tt>value</tt> is
				      * not a finite number. Still, it
				      * is allowed to store zero
				      * values in non-existent fields.
				      */
    void add (const unsigned int i,
              const unsigned int j,
	      const number value);

				     /**
				      * Multiply the entire matrix by a
				      * fixed factor.
				      */
    ChunkSparseMatrix & operator *= (const number factor);
    
				     /**
				      * Divide the entire matrix by a
				      * fixed factor.
				      */
    ChunkSparseMatrix & operator /= (const number factor);
    
				     /**
				      * Symmetrize the matrix by
				      * forming the mean value between
				      * the existing matrix and its
				      * transpose, $A = \frac 12(A+A^T)$.
				      *
				      * This operation assumes that
				      * the underlying sparsity
				      * pattern represents a symmetric
				      * object. If this is not the
				      * case, then the result of this
				      * operation will not be a
				      * symmetric matrix, since it
				      * only explicitly symmetrizes
				      * by looping over the lower left
				      * triangular part for efficiency
				      * reasons; if there are entries
				      * in the upper right triangle,
				      * then these elements are missed
				      * in the
				      * symmetrization. Symmetrization
				      * of the sparsity pattern can be
				      * obtain by
				      * ChunkSparsityPattern::symmetrize().
				      */
    void symmetrize ();
    
				     /**
				      * Copy the given matrix to this
				      * one.  The operation throws an
				      * error if the sparsity patterns
				      * of the two involved matrices
				      * do not point to the same
				      * object, since in this case the
				      * copy operation is
				      * cheaper. Since this operation
				      * is notheless not for free, we
				      * do not make it available
				      * through <tt>operator =</tt>,
				      * since this may lead to
				      * unwanted usage, e.g. in copy
				      * arguments to functions, which
				      * should really be arguments by
				      * reference.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      *
				      * The function returns a reference to
				      * <tt>*this</tt>.
				      */
    template <typename somenumber>
    ChunkSparseMatrix<number> &
    copy_from (const ChunkSparseMatrix<somenumber> &source);

				     /**
				      * This function is complete
				      * analogous to the
				      * ChunkSparsityPattern::copy_from()
				      * function in that it allows to
				      * initialize a whole matrix in
				      * one step. See there for more
				      * information on argument types
				      * and their meaning. You can
				      * also find a small example on
				      * how to use this function
				      * there.
				      *
				      * The only difference to the
				      * cited function is that the
				      * objects which the inner
				      * iterator points to need to be
				      * of type <tt>std::pair<unsigned
				      * int, value</tt>, where
				      * <tt>value</tt> needs to be
				      * convertible to the element
				      * type of this class, as
				      * specified by the
				      * <tt>number</tt> template
				      * argument.
				      *
				      * Previous content of the matrix
				      * is overwritten. Note that the
				      * entries specified by the input
				      * parameters need not
				      * necessarily cover all elements
				      * of the matrix. Elements not
				      * covered remain untouched.
				      */
    template <typename ForwardIterator>
    void copy_from (const ForwardIterator begin,
		    const ForwardIterator end);    

				     /**
				      * Copy the nonzero entries of a
				      * full matrix into this
				      * object. Previous content is
				      * deleted. Note that the
				      * underlying sparsity pattern
				      * must be appropriate to hold
				      * the nonzero entries of the
				      * full matrix.
				      */
    template <typename somenumber>
    void copy_from (const FullMatrix<somenumber> &matrix);
    
				     /**
				      * Add <tt>matrix</tt> scaled by
				      * <tt>factor</tt> to this matrix,
				      * i.e. the matrix <tt>factor*matrix</tt>
				      * is added to <tt>this</tt>. This
				      * function throws an error if the
				      * sparsity patterns of the two involved
				      * matrices do not point to the same
				      * object, since in this case the
				      * operation is cheaper.
				      *
				      * The source matrix may be a sparse
				      * matrix over an arbitrary underlying
				      * scalar type, as long as its data type
				      * is convertible to the data type of
				      * this matrix.
				      */
    template <typename somenumber>
    void add (const number factor,
	      const ChunkSparseMatrix<somenumber> &matrix);    
    
//@}
/**
 * @name Entry Access
 */
//@{
    
				     /**
				      * Return the value of the entry
				      * (<i>i,j</i>).  This may be an
				      * expensive operation and you
				      * should always take care where
				      * to call this function.  In
				      * order to avoid abuse, this
				      * function throws an exception
				      * if the required element does
				      * not exist in the matrix.
				      *
				      * In case you want a function
				      * that returns zero instead (for
				      * entries that are not in the
				      * sparsity pattern of the
				      * matrix), use the el()
				      * function.
				      *
				      * If you are looping over all elements,
				      * consider using one of the iterator
				      * classes instead, since they are
				      * tailored better to a sparse matrix
				      * structure.
				      */
    number operator () (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * This function is mostly like
				      * operator()() in that it
				      * returns the value of the
				      * matrix entry (<i>i,j</i>). The
				      * only difference is that if
				      * this entry does not exist in
				      * the sparsity pattern, then
				      * instead of raising an
				      * exception, zero is
				      * returned. While this may be
				      * convenient in some cases, note
				      * that it is simple to write
				      * algorithms that are slow
				      * compared to an optimal
				      * solution, since the sparsity
				      * of the matrix is not used.
				      *
				      * If you are looping over all elements,
				      * consider using one of the iterator
				      * classes instead, since they are
				      * tailored better to a sparse matrix
				      * structure.
				      */
    number el (const unsigned int i,
	       const unsigned int j) const;

				     /**
				      * Return the main diagonal
				      * element in the <i>i</i>th
				      * row. This function throws an
				      * error if the matrix is not
				      * quadratic (see
				      * ChunkSparsityPattern::optimize_diagonal()).
				      *
				      * This function is considerably
				      * faster than the operator()(),
				      * since for quadratic matrices, the
				      * diagonal entry may be the
				      * first to be stored in each row
				      * and access therefore does not
				      * involve searching for the
				      * right column number.
				      */
    number diag_element (const unsigned int i) const;

				     /**
				      * Same as above, but return a
				      * writeable reference. You're
				      * sure you know what you do?
				      */
    number & diag_element (const unsigned int i);

//@}
/**
 * @name Matrix vector multiplications
 */
//@{
				     /**
				      * Matrix-vector multiplication:
				      * let <i>dst = M*src</i> with
				      * <i>M</i> being this matrix.
                                      *
				      * Note that while this function can
				      * operate on all vectors that offer
				      * iterator classes, it is only really
				      * effective for objects of type @ref
				      * Vector. For all classes for which
				      * iterating over elements, or random
				      * member access is expensive, this
				      * function is not efficient. In
				      * particular, if you want to multiply
				      * with BlockVector objects, you should
				      * consider using a BlockChunkSparseMatrix as
				      * well.
				      * 
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template <class OutVector, class InVector>
    void vmult (OutVector& dst,
		const InVector& src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let <i>dst = M<sup>T</sup>*src</i> with
				      * <i>M</i> being this
				      * matrix. This function does the
				      * same as vmult() but takes
				      * the transposed matrix.
                                      *
				      * Note that while this function can
				      * operate on all vectors that offer
				      * iterator classes, it is only really
				      * effective for objects of type @ref
				      * Vector. For all classes for which
				      * iterating over elements, or random
				      * member access is expensive, this
				      * function is not efficient. In
				      * particular, if you want to multiply
				      * with BlockVector objects, you should
				      * consider using a BlockChunkSparseMatrix as
				      * well.
				      * 
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template <class OutVector, class InVector>
    void Tvmult (OutVector& dst,
		 const InVector& src) const;
  
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add
				      * <i>M*src</i> on <i>dst</i>
				      * with <i>M</i> being this
				      * matrix.
                                      *
				      * Note that while this function can
				      * operate on all vectors that offer
				      * iterator classes, it is only really
				      * effective for objects of type @ref
				      * Vector. For all classes for which
				      * iterating over elements, or random
				      * member access is expensive, this
				      * function is not efficient. In
				      * particular, if you want to multiply
				      * with BlockVector objects, you should
				      * consider using a BlockChunkSparseMatrix as
				      * well.
				      * 
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template <class OutVector, class InVector>
    void vmult_add (OutVector& dst,
		    const InVector& src) const;
    
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add
				      * <i>M<sup>T</sup>*src</i> to
				      * <i>dst</i> with <i>M</i> being
				      * this matrix. This function
				      * does the same as vmult_add()
				      * but takes the transposed
				      * matrix.
                                      *
				      * Note that while this function can
				      * operate on all vectors that offer
				      * iterator classes, it is only really
				      * effective for objects of type @ref
				      * Vector. For all classes for which
				      * iterating over elements, or random
				      * member access is expensive, this
				      * function is not efficient. In
				      * particular, if you want to multiply
				      * with BlockVector objects, you should
				      * consider using a BlockChunkSparseMatrix as
				      * well.
				      * 
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template <class OutVector, class InVector>
    void Tvmult_add (OutVector& dst,
		     const InVector& src) const;
  
				     /**
				      * Return the square of the norm
				      * of the vector $v$ with respect
				      * to the norm induced by this
				      * matrix,
				      * i.e. $\left(v,Mv\right)$. This
				      * is useful, e.g. in the finite
				      * element context, where the
				      * $L_2$ norm of a function
				      * equals the matrix norm with
				      * respect to the mass matrix of
				      * the vector representing the
				      * nodal values of the finite
				      * element function.
				      *
				      * Obviously, the matrix needs to be
				      * quadratic for this operation, and for
				      * the result to actually be a norm it
				      * also needs to be either real symmetric
				      * or complex hermitian.
				      *
				      * The underlying template types of both
				      * this matrix and the given vector
				      * should either both be real or
				      * complex-valued, but not mixed, for
				      * this function to make sense.
				      */
    template <typename somenumber>
    somenumber matrix_norm_square (const Vector<somenumber> &v) const;

				     /**
				      * Compute the matrix scalar
				      * product $\left(u,Mv\right)$.
				      */
    template <typename somenumber>
    somenumber matrix_scalar_product (const Vector<somenumber> &u,
				      const Vector<somenumber> &v) const;
				     /**
				      * Compute the residual of an
				      * equation <i>Mx=b</i>, where
				      * the residual is defined to be
				      * <i>r=b-Mx</i>. Write the
				      * residual into
				      * <tt>dst</tt>. The
				      * <i>l<sub>2</sub></i> norm of
				      * the residual vector is
				      * returned.
                                      *
                                      * Source <i>x</i> and destination
                                      * <i>dst</i> must not be the same
                                      * vector.
				      */
    template <typename somenumber>
    somenumber residual (Vector<somenumber>       &dst,
			 const Vector<somenumber> &x,
			 const Vector<somenumber> &b) const;
    
//@}
/**
 * @name Matrix norms
 */
//@{
    
    				     /**
				      * Return the l1-norm of the matrix, that is
				      * $|M|_1=max_{all columns j}\sum_{all 
				      * rows i} |M_ij|$,
				      * (max. sum of columns).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the l1-norm for vectors, i.e.
				      * $|Mv|_1\leq |M|_1 |v|_1$.
				      * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
				      */
    real_type l1_norm () const;

    				     /**
				      * Return the linfty-norm of the
				      * matrix, that is
				      * $|M|_infty=max_{all rows i}\sum_{all 
				      * columns j} |M_ij|$,
				      * (max. sum of rows).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the linfty-norm of vectors, i.e.
				      * $|Mv|_infty \leq |M|_infty |v|_infty$.
				      * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
				      */
    real_type linfty_norm () const;

                                     /**
                                      * Return the frobenius norm of the
                                      * matrix, i.e. the square root of the
                                      * sum of squares of all entries in the
                                      * matrix.
                                      */
    real_type frobenius_norm () const;
//@}
/**
 * @name Preconditioning methods
 */
//@{
    
				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the <tt>src</tt> vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * relaxation factor <tt>omega</tt>.
				      */
    template <typename somenumber>
    void precondition_Jacobi (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src,
			      const number              omega = 1.) const;

				     /**
				      * Apply SSOR preconditioning to
				      * <tt>src</tt>.
				      */
    template <typename somenumber>
    void precondition_SSOR (Vector<somenumber>       &dst,
			    const Vector<somenumber> &src,
			    const number              om = 1.) const;

				     /**
				      * Apply SOR preconditioning
				      * matrix to <tt>src</tt>.
				      */
    template <typename somenumber>
    void precondition_SOR (Vector<somenumber>       &dst,
			   const Vector<somenumber> &src,
 			   const number              om = 1.) const;
    
				     /**
				      * Apply transpose SOR
				      * preconditioning matrix to
				      * <tt>src</tt>.
				      */
    template <typename somenumber>
    void precondition_TSOR (Vector<somenumber>       &dst,
			    const Vector<somenumber> &src,
			    const number              om = 1.) const;
    
				     /**
				      * Perform SSOR preconditioning
				      * in-place.  Apply the
				      * preconditioner matrix without
				      * copying to a second vector.
				      * <tt>omega</tt> is the relaxation
				      * parameter.
				      */
    template <typename somenumber>
    void SSOR (Vector<somenumber> &v,
	       const number        omega = 1.) const;

				     /**
				      * Perform an SOR preconditioning
				      * in-place.  <tt>omega</tt> is
				      * the relaxation parameter.
				      */
    template <typename somenumber>
    void SOR (Vector<somenumber> &v,
	      const number        om = 1.) const;

				     /**
				      * Perform a transpose SOR
				      * preconditioning in-place.
				      * <tt>omega</tt> is the
				      * relaxation parameter.
				      */
    template <typename somenumber>
    void TSOR (Vector<somenumber> &v,
	       const number        om = 1.) const;

				     /**
				      * Perform a permuted SOR
				      * preconditioning in-place.
				      *
				      * The standard SOR method is
				      * applied in the order
				      * prescribed by <tt>permutation</tt>,
				      * that is, first the row
				      * <tt>permutation[0]</tt>, then
				      * <tt>permutation[1]</tt> and so
				      * on. For efficiency reasons,
				      * the permutation as well as its
				      * inverse are required.
				      *
				      * <tt>omega</tt> is the
				      * relaxation parameter.
				      */
    template <typename somenumber>
    void PSOR (Vector<somenumber> &v,
	       const std::vector<unsigned int>& permutation,
	       const std::vector<unsigned int>& inverse_permutation,
	       const number        om = 1.) const;

				     /**
				      * Perform a transposed permuted SOR
				      * preconditioning in-place.
				      *
				      * The transposed SOR method is
				      * applied in the order
				      * prescribed by
				      * <tt>permutation</tt>, that is,
				      * first the row
				      * <tt>permutation[m()-1]</tt>,
				      * then
				      * <tt>permutation[m()-2]</tt>
				      * and so on. For efficiency
				      * reasons, the permutation as
				      * well as its inverse are
				      * required.
				      *
				      * <tt>omega</tt> is the
				      * relaxation parameter.
				      */
    template <typename somenumber>
    void TPSOR (Vector<somenumber> &v,
		const std::vector<unsigned int>& permutation,
		const std::vector<unsigned int>& inverse_permutation,
		const number        om = 1.) const;

				     /**
				      * Do one SOR step on <tt>v</tt>.
				      * Performs a direct SOR step
				      * with right hand side
				      * <tt>b</tt>.
				      */
    template <typename somenumber>
    void SOR_step (Vector<somenumber> &v,
		   const Vector<somenumber> &b,
		   const number        om = 1.) const;

				     /**
				      * Do one adjoint SOR step on
				      * <tt>v</tt>.  Performs a direct
				      * TSOR step with right hand side
				      * <tt>b</tt>.
				      */
    template <typename somenumber>
    void TSOR_step (Vector<somenumber> &v,
		    const Vector<somenumber> &b,
		    const number        om = 1.) const;

				     /**
				      * Do one SSOR step on
				      * <tt>v</tt>.  Performs a direct
				      * SSOR step with right hand side
				      * <tt>b</tt> by performing TSOR
				      * after SOR.
				      */
    template <typename somenumber>
    void SSOR_step (Vector<somenumber> &v,
		    const Vector<somenumber> &b,
		    const number        om = 1.) const;
//@}
/**
 * @name Input/Output
 */
//@{

				     /**
				      * Print the matrix to the given
				      * stream, using the format
				      * <tt>(line,col) value</tt>,
				      * i.e. one nonzero entry of the
				      * matrix per line.
				      */
    void print (std::ostream &out) const;

				     /**
				      * Print the matrix in the usual
				      * format, i.e. as a matrix and
				      * not as a list of nonzero
				      * elements. For better
				      * readability, elements not in
				      * the matrix are displayed as
				      * empty space, while matrix
				      * elements which are explicitly
				      * set to zero are displayed as
				      * such.
				      *
				      * The parameters allow for a
				      * flexible setting of the output
				      * format: <tt>precision</tt> and
				      * <tt>scientific</tt> are used
				      * to determine the number
				      * format, where <tt>scientific =
				      * false</tt> means fixed point
				      * notation.  A zero entry for
				      * <tt>width</tt> makes the
				      * function compute a width, but
				      * it may be changed to a
				      * positive value, if output is
				      * crude.
				      *
				      * Additionally, a character for
				      * an empty value may be
				      * specified.
				      *
				      * Finally, the whole matrix can
				      * be multiplied with a common
				      * denominator to produce more
				      * readable output, even
				      * integers.
				      *
				      * @attention This function may
				      * produce <b>large</b> amounts
				      * of output if applied to a
				      * large matrix!
				      */
    void print_formatted (std::ostream       &out,
			  const unsigned int  precision   = 3,
			  const bool          scientific  = true,
			  const unsigned int  width       = 0,
			  const char         *zero_string = " ",
			  const double        denominator = 1.) const;

				     /**
				      * Print the actual pattern of
				      * the matrix. For each entry
				      * with an absolute value larger
				      * than threshold, a '*' is
				      * printed, a ':' for every value
				      * smaller and a '.' for every
				      * entry not allocated.
				      */
    void print_pattern(std::ostream& out,
		       const double threshold = 0.) const;
    
				     /**
				      * Write the data of this object
				      * en bloc to a file. This is
				      * done in a binary mode, so the
				      * output is neither readable by
				      * humans nor (probably) by other
				      * computers using a different
				      * operating system of number
				      * format.
				      *
				      * The purpose of this function
				      * is that you can swap out
				      * matrices and sparsity pattern
				      * if you are short of memory,
				      * want to communicate between
				      * different programs, or allow
				      * objects to be persistent
				      * across different runs of the
				      * program.
				      */
    void block_write (std::ostream &out) const;

				     /**
				      * Read data that has previously
				      * been written by block_write()
				      * from a file. This is done
				      * using the inverse operations
				      * to the above function, so it
				      * is reasonably fast because the
				      * bitstream is not interpreted
				      * except for a few numbers up
				      * front.
				      *
				      * The object is resized on this
				      * operation, and all previous
				      * contents are lost. Note,
				      * however, that no checks are
				      * performed whether new data and
				      * the underlying ChunkSparsityPattern
				      * object fit together. It is
				      * your responsibility to make
				      * sure that the sparsity pattern
				      * and the data to be read match.
				      *
				      * A primitive form of error
				      * checking is performed which
				      * will recognize the bluntest
				      * attempts to interpret some
				      * data as a matrix stored
				      * bitwise to a file that wasn't
				      * actually created that way, but
				      * not more.
				      */
    void block_read (std::istream &in);
//@}
    				     /** @addtogroup Exceptions
				      * @{ */

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The entry with index <" << arg1 << ',' << arg2
		    << "> does not exist.");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidIndex1,
		    int,
		    << "The index " << arg1 << " is not in the allowed range.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcDifferentChunkSparsityPatterns);
				     /**
				      * Exception
				      */
    DeclException2 (ExcIteratorRange,
		    int, int,
		    << "The iterators denote a range of " << arg1
		    << " elements, but the given number of rows was " << arg2);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcSourceEqualsDestination);
				     //@}    
  private:
				     /**
				      * Pointer to the sparsity
				      * pattern used for this
				      * matrix. In order to guarantee
				      * that it is not deleted while
				      * still in use, we subscribe to
				      * it using the SmartPointer
				      * class.
				      */
    SmartPointer<const ChunkSparsityPattern,ChunkSparseMatrix<number> > cols;
    
				     /**
				      * Array of values for all the
				      * nonzero entries. The position
				      * within the matrix, i.e.  the
				      * row and column number for a
				      * given entry can only be
				      * deduced using the sparsity
				      * pattern. The same holds for
				      * the more common operation of
				      * finding an entry by its
				      * coordinates.
				      */
    number *val;

				     /**
				      * Allocated size of #val. This
				      * can be larger than the
				      * actually used part if the size
				      * of the matrix was reduced
				      * somewhen in the past by
				      * associating a sparsity pattern
				      * with a smaller size to this
				      * object, using the reinit()
				      * function.
				      */
    unsigned int max_len;
    
				     /**
				      * Return the location of entry
				      * $(i,j)$ within the val array.
				      */
    unsigned int compute_location (const unsigned int i,
				   const unsigned int j) const;
    
				     // make all other sparse matrices
				     // friends
    template <typename somenumber> friend class ChunkSparseMatrix;
};

/*@}*/

#ifndef DOXYGEN
/*---------------------- Inline functions -----------------------------------*/



template <typename number>
inline
unsigned int ChunkSparseMatrix<number>::m () const
{
  Assert (cols != 0, ExcNotInitialized());
  return cols->rows;
}


template <typename number>
inline
unsigned int ChunkSparseMatrix<number>::n () const
{
  Assert (cols != 0, ExcNotInitialized());
  return cols->cols;
}



template <typename number>
inline
unsigned int
ChunkSparseMatrix<number>::compute_location (const unsigned int i,
					     const unsigned int j) const
{
  const unsigned int chunk_size = cols->get_chunk_size();
  const unsigned int chunk_index
    = cols->sparsity_pattern(i/chunk_size, j/chunk_size);

  if (chunk_index == ChunkSparsityPattern::invalid_entry)
    return ChunkSparsityPattern::invalid_entry;
  else
    {
      return (chunk_index * chunk_size * chunk_size
	      +
	      (i % chunk_size) * chunk_size
	      +
	      (j % chunk_size));
    }
}


template <typename number>
inline
void ChunkSparseMatrix<number>::set (const unsigned int i,
				     const unsigned int j,
				     const number value)
{

  Assert (numbers::is_finite(value), ExcNumberNotFinite());

  Assert (cols != 0, ExcNotInitialized());
				   // it is allowed to set elements of
				   // the matrix that are not part of
				   // the sparsity pattern, if the
				   // value to which we set it is zero
  const unsigned int index = compute_location(i,j);
  Assert ((index != SparsityPattern::invalid_entry) ||
	  (value == 0.),
	  ExcInvalidIndex(i,j));

  if (index != SparsityPattern::invalid_entry)
    val[index] = value;
}



template <typename number>
inline
void ChunkSparseMatrix<number>::add (const unsigned int i,
				     const unsigned int j,
				     const number value)
{

  Assert (numbers::is_finite(value), ExcNumberNotFinite());

  Assert (cols != 0, ExcNotInitialized());

  const unsigned int index = compute_location(i,j);
  Assert ((index != ChunkSparsityPattern::invalid_entry) ||
	  (value == 0.),
	  ExcInvalidIndex(i,j));

  if (value != 0.)
    val[index] += value;
}



template <typename number>
inline
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator *= (const number factor)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  const unsigned int chunk_size = cols->get_chunk_size();

				   // multiply all elements of the matrix with
				   // the given factor. this includes the
				   // padding elements in chunks that overlap
				   // the boundaries of the actual matrix --
				   // but since multiplication with a number
				   // does not violate the invariant of
				   // keeping these elements at zero nothing
				   // can happen
  number             *val_ptr    = val;
  const number *const end_ptr    = val +
				   cols->sparsity_pattern.n_nonzero_elements()
				   *
				   chunk_size * chunk_size;
  while (val_ptr != end_ptr)
    *val_ptr++ *= factor;

  return *this;
}



template <typename number>
inline
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator /= (const number factor)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (factor !=0, ExcDivideByZero());

  const number factor_inv = 1. / factor;

  const unsigned int chunk_size = cols->get_chunk_size();

				   // multiply all elements of the matrix with
				   // the given factor. this includes the
				   // padding elements in chunks that overlap
				   // the boundaries of the actual matrix --
				   // but since multiplication with a number
				   // does not violate the invariant of
				   // keeping these elements at zero nothing
				   // can happen
  number             *val_ptr    = val;
  const number *const end_ptr    = val +
				   cols->sparsity_pattern.n_nonzero_elements()
				   *
				   chunk_size * chunk_size;

  while (val_ptr != end_ptr)
    *val_ptr++ *= factor_inv;

  return *this;
}



template <typename number>
inline
number ChunkSparseMatrix<number>::operator () (const unsigned int i,
					       const unsigned int j) const
{
  Assert (cols != 0, ExcNotInitialized());
  AssertThrow (compute_location(i,j) != SparsityPattern::invalid_entry,
               ExcInvalidIndex(i,j));
  return val[compute_location(i,j)];
}



template <typename number>
inline
number ChunkSparseMatrix<number>::el (const unsigned int i,
				      const unsigned int j) const
{
  Assert (cols != 0, ExcNotInitialized());
  const unsigned int index = compute_location(i,j);

  if (index != ChunkSparsityPattern::invalid_entry)
    return val[index];
  else
    return 0;
}



template <typename number>
inline
number ChunkSparseMatrix<number>::diag_element (const unsigned int i) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),  ExcNotQuadratic());
  Assert (i<m(), ExcInvalidIndex1(i));
  
				   // Use that the first element in each row
				   // of a quadratic matrix is the main
				   // diagonal of the chunk sparsity pattern
  const unsigned int chunk_size = cols->get_chunk_size();
  return val[cols->sparsity_pattern.rowstart[i/chunk_size]
	     *
	     chunk_size * chunk_size
	     +
	     (i % chunk_size) * chunk_size
	     +
	     (i % chunk_size)];
}



template <typename number>
inline
number & ChunkSparseMatrix<number>::diag_element (const unsigned int i)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),  ExcNotQuadratic());
  Assert (i<m(), ExcInvalidIndex1(i));
  
				   // Use that the first element in each row
				   // of a quadratic matrix is the main
				   // diagonal of the chunk sparsity pattern
  const unsigned int chunk_size = cols->get_chunk_size();
  return val[cols->sparsity_pattern.rowstart[i/chunk_size]
	     *
	     chunk_size * chunk_size
	     +
	     (i % chunk_size) * chunk_size
	     +
	     (i % chunk_size)];
}



template <typename number>
template <typename ForwardIterator>
void
ChunkSparseMatrix<number>::copy_from (const ForwardIterator begin,
				      const ForwardIterator end)
{
  Assert (static_cast<unsigned int>(std::distance (begin, end)) == m(),
	  ExcIteratorRange (std::distance (begin, end), m()));

				   // for use in the inner loop, we
				   // define a typedef to the type of
				   // the inner iterators
  typedef typename std::iterator_traits<ForwardIterator>::value_type::const_iterator inner_iterator;
  unsigned int row=0;
  for (ForwardIterator i=begin; i!=end; ++i, ++row)
    {
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j=i->begin(); j!=end_of_row; ++j)
					 // write entries
	set (row, j->first, j->second);
    }
}


#endif // DOXYGEN


/*----------------------------   chunk_sparse_matrix.h     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   chunk_sparse_matrix.h     ---------------------------*/
