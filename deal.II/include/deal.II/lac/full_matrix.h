//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__full_matrix_h
#define __deal2__full_matrix_h


#include <deal.II/base/config.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/table.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/base/tensor.h>

#include <vector>
#include <iomanip>
#include <cstring>

DEAL_II_NAMESPACE_OPEN


// forward declarations
template <typename number> class Vector;
template <typename number> class LAPACKFullMatrix;


/*! @addtogroup Matrix1
 *@{
 */


/**
 * Implementation of a classical rectangular scheme of numbers. The
 * data type of the entries is provided in the template argument
 * <tt>number</tt>.  The interface is quite fat and in fact has grown every
 * time a new feature was needed. So, a lot of functions are provided.
 *
 * Internal calculations are usually done with the accuracy of the
 * vector argument to functions. If there is no argument with a number
 * type, the matrix number type is used.
 *
 * @note Instantiations for this template are provided for
 * <tt>@<float@>, @<double@>, @<long double@>,
 * @<std::complex@<float@>@>, @<std::complex@<double@>@>,
 * @<std::complex@<long double@>@></tt>; others can be generated in
 * application programs (see the section on @ref Instantiations in the
 * manual).
 *
 * @author Guido Kanschat, Franz-Theo Suttmeier, Wolfgang Bangerth, 1993-2004
 */
template<typename number>
class FullMatrix : public Table<2,number>
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


    class const_iterator;

				     /**
				      * Accessor class for iterators
				      */
    class Accessor
    {
      public:
					 /**
					  * Constructor. Since we use
					  * accessors only for read
					  * access, a const matrix
					  * pointer is sufficient.
					  */
	Accessor (const FullMatrix<number> *matrix,
		  const unsigned int row,
		  const unsigned int col);

					 /**
					  * Row number of the element
					  * represented by this
					  * object.
					  */
	unsigned int row() const;

					 /**
					  * Column number of the
					  * element represented by
					  * this object.
					  */
	unsigned int column() const;

					 /**
					  * Value of this matrix entry.
					  */
	number value() const;

      protected:
					 /**
					  * The matrix accessed.
					  */
	const FullMatrix<number>* matrix;

					 /**
					  * Current row number.
					  */
	unsigned int a_row;

					 /**
					  * Current column number.
					  */
	unsigned short a_col;

					 /*
					  * Make enclosing class a
					  * friend.
					  */
	friend class const_iterator;
    };

				     /**
				      * STL conforming iterator.
				      */
    class const_iterator
    {
      public:
                                         /**
                                          * Constructor.
                                          */
	const_iterator(const FullMatrix<number> *matrix,
		       const unsigned int row,
		       const unsigned int col);

                                         /**
                                          * Prefix increment.
                                          */
	const_iterator& operator++ ();

                                         /**
                                          * Postfix increment.
                                          */
	const_iterator& operator++ (int);

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor& operator* () const;

                                         /**
                                          * Dereferencing operator.
                                          */
	const Accessor* operator-> () const;

                                         /**
                                          * Comparison. True, if
                                          * both iterators point to
                                          * the same matrix
                                          * position.
                                          */
	bool operator == (const const_iterator&) const;
                                         /**
                                          * Inverse of <tt>==</tt>.
                                          */
	bool operator != (const const_iterator&) const;

                                         /**
                                          * Comparison operator. Result is
                                          * true if either the first row
                                          * number is smaller or if the row
                                          * numbers are equal and the first
                                          * index is smaller.
                                          */
	bool operator < (const const_iterator&) const;

                                         /**
                                          * Comparison operator. Compares just
                                          * the other way around than the
                                          * operator above.
                                          */
	bool operator > (const const_iterator&) const;

      private:
                                         /**
                                          * Store an object of the
                                          * accessor class.
                                          */
        Accessor accessor;
    };
/**
 * @name Constructors and initalization.
 * See also the base class Table.
 */
//@{

				     /**
				      * Constructor. Initialize the
				      * matrix as a square matrix with
				      * dimension <tt>n</tt>.
				      *
				      * In order to avoid the implicit
				      * conversion of integers and
				      * other types to a matrix, this
				      * constructor is declared
				      * <tt>explicit</tt>.
				      *
				      * By default, no memory is
				      * allocated.
				      */
    explicit FullMatrix (const unsigned int n = 0);

				     /**
				      * Constructor. Initialize the
				      * matrix as a rectangular
				      * matrix.
				      */
    FullMatrix (const unsigned int rows,
                const unsigned int cols);

				     /**
				      * Copy constructor. This
				      * constructor does a deep copy
				      * of the matrix. Therefore, it
				      * poses a possible efficiency
				      * problem, if for example,
				      * function arguments are passed
				      * by value rather than by
				      * reference. Unfortunately, we
				      * can't mark this copy
				      * constructor <tt>explicit</tt>,
				      * since that prevents the use of
				      * this class in containers, such
				      * as <tt>std::vector</tt>. The
				      * responsibility to check
				      * performance of programs must
				      * therefore remain with the
				      * user of this class.
				      */
    FullMatrix (const FullMatrix&);

				     /**
				      * Constructor initializing from
				      * an array of numbers. The array
				      * is arranged line by line. No
				      * range checking is performed.
				      */
    FullMatrix (const unsigned int rows,
		const unsigned int cols,
		const number* entries);

				     /**
				      * Construct a full matrix that
				      * equals the identity matrix of
				      * the size of the
				      * argument. Using this
				      * constructor, one can easily
				      * create an identity matrix of
				      * size <code>n</code> by saying
				      * @verbatim
				      * FullMatrix<double> M(IdentityMatrix(n));
				      * @endverbatim
				      */
    explicit FullMatrix (const IdentityMatrix &id);


				     /**
				      * Assignment operator.
				      */
    FullMatrix<number> &
    operator = (const FullMatrix<number>&);

				     /**
				      * Variable assignment operator.
				      */
    template<typename number2>
    FullMatrix<number> &
    operator = (const FullMatrix<number2>&);

				     /**
				      * This operator assigns a scalar
				      * to a matrix. To avoid
				      * confusion with the semantics
				      * of this function, zero is the
				      * only value allowed for
				      * <tt>d</tt>, allowing you to
				      * clear a matrix in an intuitive
				      * way.
				      */
    FullMatrix<number> &
    operator = (const number d);

				     /**
				      * Copy operator to create a full
				      * matrix that equals the
				      * identity matrix of the size of
				      * the argument. This way, one can easily
				      * create an identity matrix of
				      * size <code>n</code> by saying
				      * @verbatim
				      *   M = IdentityMatrix(n);
				      * @endverbatim
				      */
    FullMatrix<number> &
    operator = (const IdentityMatrix &id);

				     /**
				      * Assignment operator for a
				      * LapackFullMatrix. The calling matrix
				      * must be of the same size as the
				      * LAPACK matrix.
				      */
    template <typename number2>
    FullMatrix<number> &
    operator = (const LAPACKFullMatrix<number2>&);

                                     /**
				      * Assignment from different
				      * matrix classes. This
				      * assignment operator uses
				      * iterators of the class
				      * MATRIX. Therefore, sparse
				      * matrices are possible sources.
				      */
    template <class MATRIX>
    void copy_from (const MATRIX&);

				     /**
				      * Fill rectangular block.
				      *
				      * A rectangular block of the
				      * matrix <tt>src</tt> is copied into
				      * <tt>this</tt>. The upper left
				      * corner of the block being
				      * copied is
				      * <tt>(src_offset_i,src_offset_j)</tt>.
				      * The upper left corner of the
				      * copied block is
				      * <tt>(dst_offset_i,dst_offset_j)</tt>.
				      * The size of the rectangular
				      * block being copied is the
				      * maximum size possible,
				      * determined either by the size
				      * of <tt>this</tt> or <tt>src</tt>.
				      */
    template<typename number2>
    void fill (const FullMatrix<number2> &src,
	       const unsigned int dst_offset_i = 0,
	       const unsigned int dst_offset_j = 0,
	       const unsigned int src_offset_i = 0,
	       const unsigned int src_offset_j = 0);


				     /**
				      * Make function of base class
				      * available.
				      */
    template<typename number2>
    void fill (const number2*);

				     /**
				      * Fill with permutation of
				      * another matrix.
				      *
				      * The matrix <tt>src</tt> is copied
				      * into the target. The two
				      * permutation <tt>p_r</tt> and
				      * <tt>p_c</tt> operate in a way, such
				      * that <tt>result(i,j) =
				      * src(p_r[i], p_c[j])</tt>.
				      *
				      * The vectors may also be a
				      * selection from a larger set of
				      * integers, if the matrix
				      * <tt>src</tt> is bigger. It is also
				      * possible to duplicate rows or
				      * columns by this method.
				      */
    template<typename number2>
    void fill_permutation (const FullMatrix<number2>       &src,
			   const std::vector<unsigned int> &p_rows,
			   const std::vector<unsigned int> &p_cols);

//@}
///@name Non-modifying operators
//@{

				     /**
				      * Comparison operator. Be
				      * careful with this thing, it
				      * may eat up huge amounts of
				      * computing time! It is most
				      * commonly used for internal
				      * consistency checks of
				      * programs.
				      */
    bool operator == (const FullMatrix<number> &) const;

				     /**
				      * Number of rows of this matrix.
				      * To remember: this matrix is an
				      * <i>m x n</i>-matrix.
				      */
    unsigned int m () const;

				     /**
				      * Number of columns of this matrix.
				      * To remember: this matrix is an
				      * <i>m x n</i>-matrix.
				      */
    unsigned int n () const;

    				     /**
				      * Return whether the matrix
				      * contains only elements with
				      * value zero. This function is
				      * mainly for internal
				      * consistency checks and should
				      * seldomly be used when not in
				      * debug mode since it uses quite
				      * some time.
				      */
    bool all_zero () const;

    				     /**
				      * Return the square of the norm
				      * of the vector <tt>v</tt> induced by
				      * this matrix,
				      * i.e. <i>(v,Mv)</i>. This is
				      * useful, e.g. in the finite
				      * element context, where the
				      * <i>L<sup>2</sup></i> norm of a
				      * function equals the matrix
				      * norm with respect to the mass
				      * matrix of the vector
				      * representing the nodal values
				      * of the finite element
				      * function.
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
    template<typename number2>
    number2 matrix_norm_square (const Vector<number2> &v) const;

				     /**
				      * Build the matrix scalar
				      * product <tt>u<sup>T</sup> M
				      * v</tt>. This function is
				      * mostly useful when building
				      * the cellwise scalar product of
				      * two functions in the finite
				      * element context.
				      *
				      * The underlying template types of both
				      * this matrix and the given vector
				      * should either both be real or
				      * complex-valued, but not mixed, for
				      * this function to make sense.
				      */
    template<typename number2>
    number2 matrix_scalar_product (const Vector<number2> &u,
				   const Vector<number2> &v) const;

    				     /**
				      * Return the
				      * <i>l<sub>1</sub></i>-norm of
				      * the matrix, where
				      * $||M||_1 = \max_j \sum_i
				      * |M_{ij}|$ (maximum of
				      * the sums over columns).
				      */
    real_type l1_norm () const;

    				     /**
				      * Return the $l_\infty$-norm of
				      * the matrix, where
				      * $||M||_\infty = \max_i \sum_j
				      * |M_{ij}|$ (maximum of the sums
				      * over rows).
				      */
    real_type linfty_norm () const;

				     /**
				      * Compute the Frobenius norm of
				      * the matrix.  Return value is
				      * the root of the square sum of
				      * all matrix entries.
				      *
				      * @note For the timid among us:
				      * this norm is not the norm
				      * compatible with the
				      * <i>l<sub>2</sub></i>-norm of
				      * the vector space.
				      */
    real_type frobenius_norm () const;

				     /**
				      * Compute the relative norm of
				      * the skew-symmetric part. The
				      * return value is the Frobenius
				      * norm of the skew-symmetric
				      * part of the matrix divided by
				      * that of the matrix.
				      *
				      * Main purpose of this function
				      * is to check, if a matrix is
				      * symmetric within a certain
				      * accuracy, or not.
				      */
    real_type relative_symmetry_norm2 () const;

				     /**
                                      * Computes the determinant of a
                                      * matrix.  This is only
                                      * implemented for one, two, and
                                      * three dimensions, since for
                                      * higher dimensions the
                                      * numerical work explodes.
                                      * Obviously, the matrix needs to
                                      * be quadratic for this function.
                                      */
    number determinant () const;

				     /**
				      * Return the trace of the matrix,
				      * i.e. the sum of the diagonal values
				      * (which happens to also equal the sum
				      * of the eigenvalues of a matrix).
                                      * Obviously, the matrix needs to
                                      * be quadratic for this function.
                                      */
    number trace () const;

				     /**
				      * Output of the matrix in
				      * user-defined format.
				      */
    template <class STREAM>
    void print (STREAM             &s,
		const unsigned int  width=5,
		const unsigned int  precision=2) const;

				     /**
				      * Print the matrix and allow
				      * formatting of entries.
				      *
				      * The parameters allow for a
				      * flexible setting of the output
				      * format:
				      *
				      * @arg <tt>precision</tt>
				      * denotes the number of trailing
				      * digits.
				      *
				      * @arg <tt>scientific</tt> is
				      * used to determine the number
				      * format, where
				      * <tt>scientific</tt> =
				      * <tt>false</tt> means fixed
				      * point notation.
				      *
				      * @arg <tt>width</tt> denotes
				      * the with of each column. A
				      * zero entry for <tt>width</tt>
				      * makes the function compute a
				      * width, but it may be changed
				      * to a positive value, if output
				      * is crude.
				      *
				      * @arg <tt>zero_string</tt>
				      * specifies a string printed for
				      * zero entries.
				      *
				      * @arg <tt>denominator</tt>
				      * Multiply the whole matrix by
				      * this common denominator to get
				      * nicer numbers.
				      *
				      * @arg <tt>threshold</tt>: all
				      * entries with absolute value
				      * smaller than this are
				      * considered zero.
				     */
    void print_formatted (std::ostream       &out,
			  const unsigned int  precision=3,
			  const bool          scientific  = true,
			  const unsigned int  width       = 0,
			  const char         *zero_string = " ",
			  const double        denominator = 1.,
			  const double        threshold   = 0.) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    std::size_t memory_consumption () const;

//@}
///@name Iterator functions
//@{

				     /**
				      * STL-like iterator with the
				      * first entry.
				      */
    const_iterator begin () const;

				     /**
				      * Final iterator.
				      */
    const_iterator end () const;

				     /**
				      * STL-like iterator with the
				      * first entry of row <tt>r</tt>.
				      */
    const_iterator begin (const unsigned int r) const;

				     /**
				      * Final iterator of row <tt>r</tt>.
				      */
    const_iterator end (const unsigned int r) const;

//@}
///@name Modifying operators
//@{

				     /**
				      * Scale the entire matrix by a
				      * fixed factor.
				      */
    FullMatrix & operator *= (const number factor);

				     /**
				      * Scale the entire matrix by the
				      * inverse of the given factor.
				      */
    FullMatrix & operator /= (const number factor);

                                     /**
				      * Simple addition of a scaled
				      * matrix, i.e. <tt>*this +=
				      * a*A</tt>.
				      *
				      * The matrix <tt>A</tt> may be a
				      * full matrix over an arbitrary
				      * underlying scalar type, as
				      * long as its data type is
				      * convertible to the data type
				      * of this matrix.
				      */
    template<typename number2>
    void add (const number               a,
	      const FullMatrix<number2> &A);

                                     /**
				      * Multiple addition of scaled
				      * matrices, i.e. <tt>*this +=
				      * a*A + b*B</tt>.
				      *
				      * The matrices <tt>A</tt> and
				      * <tt>B</tt> may be a full
				      * matrix over an arbitrary
				      * underlying scalar type, as
				      * long as its data type is
				      * convertible to the data type
				      * of this matrix.
				      */
    template<typename number2>
    void add (const number               a,
	      const FullMatrix<number2> &A,
	      const number               b,
	      const FullMatrix<number2> &B);

                                     /**
				      * Multiple addition of scaled
				      * matrices, i.e. <tt>*this +=
				      * a*A + b*B + c*C</tt>.
				      *
				      * The matrices <tt>A</tt>,
				      * <tt>B</tt> and <tt>C</tt> may
				      * be a full matrix over an
				      * arbitrary underlying scalar
				      * type, as long as its data type
				      * is convertible to the data
				      * type of this matrix.
				      */
    template<typename number2>
    void add (const number               a,
	      const FullMatrix<number2> &A,
	      const number               b,
	      const FullMatrix<number2> &B,
	      const number               c,
	      const FullMatrix<number2> &C);

				     /**
				      * Add rectangular block.
				      *
				      * A rectangular block of the matrix
				      * <tt>src</tt> is added to
				      * <tt>this</tt>. The upper left corner
				      * of the block being copied is
				      * <tt>(src_offset_i,src_offset_j)</tt>.
				      * The upper left corner of the copied
				      * block is
				      * <tt>(dst_offset_i,dst_offset_j)</tt>.
				      * The size of the rectangular block
				      * being copied is the maximum size
				      * possible, determined either by the
				      * size of <tt>this</tt> or <tt>src</tt>
				      * and the given offsets.
				      */
    template<typename number2>
    void add (const FullMatrix<number2> &src,
	      const number factor,
	      const unsigned int dst_offset_i = 0,
	      const unsigned int dst_offset_j = 0,
	      const unsigned int src_offset_i = 0,
	      const unsigned int src_offset_j = 0);

				     /**
				      * Weighted addition of the
				      * transpose of <tt>B</tt> to
				      * <tt>this</tt>.
				      *
				      * <i>A += s B<sup>T</sup></i>
				      */
    template<typename number2>
    void Tadd (const number               s,
	       const FullMatrix<number2> &B);

				     /**
				      * Add transpose of a rectangular block.
				      *
				      * A rectangular block of the
				      * matrix <tt>src</tt> is
				      * transposed and addedadded to
				      * <tt>this</tt>. The upper left
				      * corner of the block being
				      * copied is
				      * <tt>(src_offset_i,src_offset_j)</tt>
				      * in the coordinates of the
				      * <b>non</b>-transposed matrix.
				      * The upper left corner of the
				      * copied block is
				      * <tt>(dst_offset_i,dst_offset_j)</tt>.
				      * The size of the rectangular
				      * block being copied is the
				      * maximum size possible,
				      * determined either by the size
				      * of <tt>this</tt> or
				      * <tt>src</tt>.
				      */
    template<typename number2>
    void Tadd (const FullMatrix<number2> &src,
	       const number               factor,
	       const unsigned int dst_offset_i = 0,
	       const unsigned int dst_offset_j = 0,
	       const unsigned int src_offset_i = 0,
	       const unsigned int src_offset_j = 0);

				     /**
				      * <i>A(i,1...n) +=
				      * s*A(j,1...n)</i>.  Simple
				      * addition of rows of this
				      */
    void add_row (const unsigned int i,
		  const number       s,
		  const unsigned int j);

				     /**
				      * <i>A(i,1...n) += s*A(j,1...n)
				      * + t*A(k,1...n)</i>.  Multiple
				      * addition of rows of this.
				      */
    void add_row (const unsigned int i,
		  const number s, const unsigned int j,
		  const number t, const unsigned int k);

				     /**
				      * <i>A(1...n,i) += s*A(1...n,j)</i>.
				      *  Simple addition of columns of this.
				      */
    void add_col (const unsigned int i,
		  const number       s,
		  const unsigned int j);

				     /**
				      * <i>A(1...n,i) += s*A(1...n,j)
				      * + t*A(1...n,k)</i>.  Multiple
				      * addition of columns of this.
				      */
    void add_col (const unsigned int i,
		  const number s, const unsigned int j,
		  const number t, const unsigned int k);

				     /**
				      * Swap <i>A(i,1...n) <->
				      * A(j,1...n)</i>.  Swap rows i
				      * and j of this
				      */
    void swap_row (const unsigned int i,
		   const unsigned int j);

				     /**
				      *  Swap <i>A(1...n,i) <->
				      *  A(1...n,j)</i>.  Swap columns
				      *  i and j of this
				      */
    void swap_col (const unsigned int i,
		   const unsigned int j);

				     /**
				      * Add constant to diagonal
				      * elements of this, i.e. add a
				      * multiple of the identity
				      * matrix.
				      */
    void diagadd (const number s);

				     /**
				      * Assignment <tt>*this =
				      * a*A</tt>.
				      */
    template<typename number2>
    void equ (const number               a,
	      const FullMatrix<number2> &A);

				     /**
				      * Assignment <tt>*this = a*A +
				      * b*B</tt>.
				      */
    template<typename number2>
    void equ (const number               a,
	      const FullMatrix<number2> &A,
	      const number               b,
	      const FullMatrix<number2> &B);

				     /**
				      * Assignment <tt>*this = a*A +
				      * b*B + c*C</tt>.
				      */
    template<typename number2>
    void equ (const number               a,
	      const FullMatrix<number2> &A,
	      const number               b,
	      const FullMatrix<number2> &B,
	      const number               c,
	      const FullMatrix<number2> &C);

				     /**
				      * Symmetrize the matrix by
				      * forming the mean value between
				      * the existing matrix and its
				      * transpose, <i>A =
				      * 1/2(A+A<sup>T</sup>)</i>.
				      *
				      * Obviously the matrix must be
				      * quadratic for this operation.
				      */
    void symmetrize ();

    				     /**
				      * A=Inverse(A). A must be a square matrix.
                                      * Inversion of
				      * this matrix by Gauss-Jordan
				      * algorithm with partial
				      * pivoting.  This process is
				      * well-behaved for positive
				      * definite matrices, but be
				      * aware of round-off errors in
				      * the indefinite case.
				      *
				      * In case deal.II was configured with
				      * LAPACK, the functions Xgetrf and
				      * Xgetri build an LU factorization and
				      * invert the matrix upon that
				      * factorization, providing best
				      * performance up to matrices with a
				      * few hundreds rows and columns.
				      *
				      * The numerical effort to invert
				      * an <tt>n x n</tt> matrix is of the
				      * order <tt>n**3</tt>.
				      */
    void gauss_jordan ();

				     /**
				      * Assign the inverse of the given matrix
				      * to <tt>*this</tt>. This function is
				      * hardcoded for quadratic matrices of
				      * dimension one to four. However, since
				      * the amount of code needed grows
				      * quickly, the method gauss_jordan() is
				      * invoked implicitly if the dimension is
				      * larger.
				      */
    template <typename number2>
    void invert (const FullMatrix<number2> &M);

				     /**
				      * Assign the Cholesky decomposition
				      * of the given matrix to <tt>*this</tt>.
				      * The given matrix must be symmetric
				      * positive definite.
				      *
				      * ExcMatrixNotPositiveDefinite
				      * will be thrown in the case that the
				      * matrix is not positive definite.
				      */
    template <typename number2>
    void cholesky (const FullMatrix<number2> &A);

				     /**
				      * <tt>*this(i,j)</tt> = $V(i) W(j)$
				      * where $V,W$
				      * are vectors of the same length.
				      */
    template <typename number2>
    void outer_product (const Vector<number2> &V,
			const Vector<number2> &W);

				     /**
				      * Assign the left_inverse of the given matrix
				      * to <tt>*this</tt>. The calculation being
				      * performed is <i>(A<sup>T</sup>*A)<sup>-1</sup>
				      * *A<sup>T</sup></i>.
				      */
    template <typename number2>
    void left_invert (const FullMatrix<number2> &M);

				     /**
				      * Assign the right_inverse of the given matrix
				      * to <tt>*this</tt>. The calculation being
				      * performed is <i>A<sup>T</sup>*(A*A<sup>T</sup>)
				      * <sup>-1</sup></i>.
				      */
    template <typename number2>
    void right_invert (const FullMatrix<number2> &M);

                                     /**
				      * Fill matrix with elements
				      * extracted from a tensor,
				      * taking rows included between
				      * <tt>r_i</tt> and <tt>r_j</tt>
				      * and columns between
				      * <tt>c_i</tt> and
				      * <tt>c_j</tt>. The resulting
				      * matrix is then inserted in the
				      * destination matrix at position
				      * <tt>(dst_r, dst_c)</tt> Checks
				      * on the indices are made.
				      */
    template <int dim>
    void
    copy_from (Tensor<2,dim> &T,
	       const unsigned int src_r_i=0,
	       const unsigned int src_r_j=dim-1,
	       const unsigned int src_c_i=0,
	       const unsigned int src_c_j=dim-1,
	       const unsigned int dst_r=0,
	       const unsigned int dst_c=0);

                                     /**
				      * Insert a submatrix (also
				      * rectangular) into a tensor,
				      * putting its upper left element
				      * at the specified position
				      * <tt>(dst_r, dst_c)</tt> and
				      * the other elements
				      * consequently. Default values
				      * are chosen so that no
				      * parameter needs to be specified
				      * if the size of the tensor and
				      * that of the matrix coincide.
                                      */
    template <int dim>
    void
    copy_to(Tensor<2,dim> &T,
	    const unsigned int src_r_i=0,
	    const unsigned int src_r_j=dim-1,
	    const unsigned int src_c_i=0,
	    const unsigned int src_c_j=dim-1,
	    const unsigned int dst_r=0,
	    const unsigned int dst_c=0);


//@}
///@name Multiplications
//@{

				     /**
				      * Matrix-matrix-multiplication.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>C</tt> or added
				      * to <tt>C</tt>.
				      *
				      * if (adding)
				      *  <i>C += A*B</i>
				      *
				      * if (!adding)
				      *  <i>C = A*B</i>
				      *
				      * Assumes that <tt>A</tt> and
				      * <tt>B</tt> have compatible sizes and
				      * that <tt>C</tt> already has the
				      * right size.
				      *
				      * This function uses the BLAS function
				      * Xgemm if the calling matrix has more
				      * than 15 rows and BLAS was detected
				      * during configuration. Using BLAS
				      * usually results in considerable
				      * performance gains.
				      */
    template<typename number2>
    void mmult (FullMatrix<number2>       &C,
		const FullMatrix<number2> &B,
		const bool                 adding=false) const;

				     /**
				      * Matrix-matrix-multiplication using
				      * transpose of <tt>this</tt>.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>C</tt> or added
				      * to <tt>C</tt>.
				      *
				      * if (adding)
				      *  <i>C += A<sup>T</sup>*B</i>
				      *
				      * if (!adding)
				      *  <i>C = A<sup>T</sup>*B</i>
				      *
				      * Assumes that <tt>A</tt> and
				      * <tt>B</tt> have compatible
				      * sizes and that <tt>C</tt>
				      * already has the right size.
				      *
				      * This function uses the BLAS function
				      * Xgemm if the calling matrix has more
				      * than 15 columns and BLAS was
				      * detected during configuration. Using
				      * BLAS usually results in considerable
				      * performance gains.
				      */
    template<typename number2>
    void Tmmult (FullMatrix<number2>       &C,
		 const FullMatrix<number2> &B,
		 const bool                 adding=false) const;

				     /**
				      * Matrix-matrix-multiplication using
				      * transpose of <tt>B</tt>.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>C</tt> or added
				      * to <tt>C</tt>.
				      *
				      * if (adding)
				      *  <i>C += A*B<sup>T</sup></i>
				      *
				      * if (!adding)
				      *  <i>C = A*B<sup>T</sup></i>
				      *
				      * Assumes that <tt>A</tt> and
				      * <tt>B</tt> have compatible sizes and
				      * that <tt>C</tt> already has the
				      * right size.
				      *
				      * This function uses the BLAS function
				      * Xgemm if the calling matrix has more
				      * than 15 rows and BLAS was detected
				      * during configuration. Using BLAS
				      * usually results in considerable
				      * performance gains.
				      */
    template<typename number2>
    void mTmult (FullMatrix<number2>       &C,
		 const FullMatrix<number2> &B,
		 const bool                 adding=false) const;

				     /**
				      * Matrix-matrix-multiplication using
				      * transpose of <tt>this</tt> and
				      * <tt>B</tt>.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>C</tt> or added
				      * to <tt>C</tt>.
				      *
				      * if (adding)
				      *  <i>C += A<sup>T</sup>*B<sup>T</sup></i>
				      *
				      * if (!adding)
				      *  <i>C = A<sup>T</sup>*B<sup>T</sup></i>
				      *
				      * Assumes that <tt>A</tt> and
				      * <tt>B</tt> have compatible
				      * sizes and that <tt>C</tt>
				      * already has the right size.
				      *
				      * This function uses the BLAS function
				      * Xgemm if the calling matrix has more
				      * than 15 columns and BLAS was
				      * detected during configuration. Using
				      * BLAS usually results in considerable
				      * performance gains.
				      */
    template<typename number2>
    void TmTmult (FullMatrix<number2>       &C,
		  const FullMatrix<number2> &B,
		  const bool                 adding=false) const;
    
				     /**
				      * Add to the current matrix the
				      * Schur complement <b>B
				      * A<sup>-1</sup>
				      * D</b>. Optionally, use the
				      * transposes of the matrices
				      * <b>B</b> and <b>D</b>. Note
				      * that the argument for <b>A</b>
				      * is already the inverse.
				      */
    void schur_complement(const FullMatrix<number>& Ainverse,
			  const FullMatrix<number>& B,
			  const FullMatrix<number>& D,
			  const bool transpose_B = false,
			  const bool transpose_D = false);
    
				     /**
				      * Matrix-vector-multiplication.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>w</tt> or added
				      * to <tt>w</tt>.
				      *
				      * if (adding)
				      *  <i>w += A*v</i>
				      *
				      * if (!adding)
				      *  <i>w = A*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template<typename number2>
    void vmult (Vector<number2>       &w,
		const Vector<number2> &v,
		const bool             adding=false) const;

				     /**
				      * Adding Matrix-vector-multiplication.
				      *  <i>w += A*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template<typename number2>
    void vmult_add (Vector<number2>       &w,
		    const Vector<number2> &v) const;

				     /**
				      * Transpose
				      * matrix-vector-multiplication.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>w</tt> or added
				      * to <tt>w</tt>.
				      *
				      * if (adding)
				      *  <i>w += A<sup>T</sup>*v</i>
				      *
				      * if (!adding)
				      *  <i>w = A<sup>T</sup>*v</i>
                                      *
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template<typename number2>
    void Tvmult (Vector<number2>       &w,
		 const Vector<number2> &v,
		 const bool             adding=false) const;

				     /**
				      * Adding transpose
				      * matrix-vector-multiplication.
				      *  <i>w += A<sup>T</sup>*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    template<typename number2>
    void Tvmult_add (Vector<number2>       &w,
		     const Vector<number2> &v) const;

				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the <tt>src</tt> vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * damping factor <tt>omega</tt>.
				      */
    template <typename somenumber>
    void precondition_Jacobi (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src,
			      const number              omega = 1.) const;

				     /**
				      * <i>dst=b-A*x</i>. Residual calculation,
				      * returns the <i>l<sub>2</sub></i>-norm
				      * |<i>dst</i>|.
                                      *
                                      * Source <i>x</i> and destination
                                      * <i>dst</i> must not be the same
                                      * vector.
				      */
    template<typename number2, typename number3>
    number residual (Vector<number2>       &dst,
		     const Vector<number2> &x,
		     const Vector<number3> &b) const;

				     /**
				      * Forward elimination of lower
				      * triangle.  Inverts the lower
				      * triangle of a rectangular matrix
				      * for a given right hand side.
				      *
				      * If the matrix has more columns
				      * than rows, this function only
				      * operates on the left quadratic
				      * submatrix. If there are more
				      * rows, the upper quadratic part
				      * of the matrix is considered.
				      *
				      * @note It is safe to use the
				      * same object for @p dst and @p
				      * src.
				      */
    template<typename number2>
    void forward (Vector<number2>       &dst,
		  const Vector<number2> &src) const;

				     /**
				      * Backward elimination of upper
				      * triangle.
				      *
				      * See forward()
				      *
				      * @note It is safe to use the
				      * same object for @p dst and @p
				      * src.
				      */
    template<typename number2>
    void backward (Vector<number2>       &dst,
		   const Vector<number2> &src) const;

                                     //@}

				     /** @addtogroup Exceptions
				      * @{ */

				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyMatrix);

				     /**
				      * Exception
				      */
    DeclException1 (ExcNotRegular,
		    number,
		    << "The maximal pivot is " << arg1
		    << ", which is below the threshold. The matrix may be singular.");
				     /**
				      * Exception
				      */
    DeclException3 (ExcInvalidDestination,
		    int, int, int,
		    << "Target region not in matrix: size in this direction="
		    << arg1 << ", size of new matrix=" << arg2
		    << ", offset=" << arg3);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcSourceEqualsDestination);
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotPositiveDefinite);
				     //@}

    friend class Accessor;
};

/**@}*/

#ifndef DOXYGEN
/*-------------------------Inline functions -------------------------------*/




template <typename number>
inline
unsigned int
FullMatrix<number>::m() const
{
  return this->n_rows();
}



template <typename number>
inline
unsigned int
FullMatrix<number>::n() const
{
  return this->n_cols();
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator = (const number d)
{
  Assert (d==number(0), ExcScalarAssignmentOnlyForZeroValue());

  if (this->n_elements() != 0)
    memset (&this->values[0], 0, this->n_elements()*sizeof(number));

  return *this;
}



template <typename number>
template <typename number2>
inline
void FullMatrix<number>::fill (const number2* src)
{
  Table<2,number>::fill(src);
}



template <typename number>
template <class MATRIX>
void
FullMatrix<number>::copy_from (const MATRIX& M)
{
  this->reinit (M.m(), M.n());
  const typename MATRIX::const_iterator end = M.end();
  for (typename MATRIX::const_iterator entry = M.begin();
       entry != end; ++entry)
    this->el(entry->row(), entry->column()) = entry->value();
}


template <typename number>
template<typename number2>
void
FullMatrix<number>::vmult_add (Vector<number2>       &w,
			       const Vector<number2> &v) const
{
  vmult(w, v, true);
}


template <typename number>
template<typename number2>
void
FullMatrix<number>::Tvmult_add (Vector<number2>       &w,
				const Vector<number2> &v) const
{
  Tvmult(w, v, true);
}


//---------------------------------------------------------------------------


template <typename number>
inline
FullMatrix<number>::Accessor::
Accessor (const FullMatrix<number>* matrix,
          const unsigned int r,
          const unsigned int c)
		:
		matrix(matrix),
		a_row(r),
		a_col(c)
{}


template <typename number>
inline
unsigned int
FullMatrix<number>::Accessor::row() const
{
  return a_row;
}


template <typename number>
inline
unsigned int
FullMatrix<number>::Accessor::column() const
{
  return a_col;
}


template <typename number>
inline
number
FullMatrix<number>::Accessor::value() const
{
  Assert (numbers::is_finite( matrix->el(a_row, a_col) ), ExcNumberNotFinite());
  return matrix->el(a_row, a_col);
}


template <typename number>
inline
FullMatrix<number>::const_iterator::
const_iterator(const FullMatrix<number> *matrix,
               const unsigned int r,
               const unsigned int c)
		:
		accessor(matrix, r, c)
{}


template <typename number>
inline
typename FullMatrix<number>::const_iterator &
FullMatrix<number>::const_iterator::operator++ ()
{
  Assert (accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

  ++accessor.a_col;
  if (accessor.a_col >= accessor.matrix->n())
    {
      accessor.a_col = 0;
      accessor.a_row++;
    }
  return *this;
}


template <typename number>
inline
const typename FullMatrix<number>::Accessor &
FullMatrix<number>::const_iterator::operator* () const
{
  return accessor;
}


template <typename number>
inline
const typename FullMatrix<number>::Accessor *
FullMatrix<number>::const_iterator::operator-> () const
{
  return &accessor;
}


template <typename number>
inline
bool
FullMatrix<number>::const_iterator::
operator == (const const_iterator& other) const
{
  return (accessor.row() == other.accessor.row() &&
          accessor.column() == other.accessor.column());
}


template <typename number>
inline
bool
FullMatrix<number>::const_iterator::
operator != (const const_iterator& other) const
{
  return ! (*this == other);
}


template <typename number>
inline
bool
FullMatrix<number>::const_iterator::
operator < (const const_iterator& other) const
{
  return (accessor.row() < other.accessor.row() ||
	  (accessor.row() == other.accessor.row() &&
           accessor.column() < other.accessor.column()));
}


template <typename number>
inline
bool
FullMatrix<number>::const_iterator::
operator > (const const_iterator& other) const
{
  return (other < *this);
}


template <typename number>
inline
typename FullMatrix<number>::const_iterator
FullMatrix<number>::begin () const
{
  return const_iterator(this, 0, 0);
}


template <typename number>
inline
typename FullMatrix<number>::const_iterator
FullMatrix<number>::end () const
{
  return const_iterator(this, m(), 0);
}


template <typename number>
inline
typename FullMatrix<number>::const_iterator
FullMatrix<number>::begin (const unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r, 0);
}



template <typename number>
inline
typename FullMatrix<number>::const_iterator
FullMatrix<number>::end (const unsigned int r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r+1, 0);
}



template <typename number>
template <class STREAM>
inline
void
FullMatrix<number>::print (STREAM             &s,
			   const unsigned int  w,
			   const unsigned int  p) const
{
  Assert (!this->empty(), ExcEmptyMatrix());

  for (unsigned int i=0; i<this->m(); ++i)
    {
      for (unsigned int j=0; j<this->n(); ++j)
	s << std::setw(w) << std::setprecision(p) << this->el(i,j);
      s << std::endl;
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

