//----------------------------  full_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix.h  ---------------------------
#ifndef __deal2__full_matrix_h
#define __deal2__full_matrix_h


/*----------------------------   fullmatrix.h     ---------------------------*/

#include <base/exceptions.h>
#include <base/vector2d.h>

#include <vector>

// forward declarations
template<typename number> class Vector;



/**
 * Rectangular/quadratic full matrix.
 *
 * Implementation of a classical rectangular scheme of numbers. The
 * data type of the entries is provided in the template argument
 * @p{number}.  The interface is quite fat and in fact has grown every
 * time a new feature was needed. So, a lot of functions are provided.
 *
 * Since the instantiation of this template is quite an effort,
 * standard versions are precompiled into the library. These include
 * all combinations of @p{float} and @p{double} for matrices and
 * vectors. If you need more data types, the implementation of
 * non-inline functions is in @p{fullmatrix.templates.h}. Driver files
 * are in the source tree.
 *
 * Internal calculations are usually done with the accuracy of the
 * vector argument to functions. If there is no argument with a number
 * type, the matrix number type is used.
 *
 *
 * @sect2{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Guido Kanschat, Franz-Theo Suttmeier, Wolfgang Bangerth, 1993-2001
 */
template<typename number>
class FullMatrix : public vector2d<number>
{
  public:
				     /**
				      * Constructor. Initialize the matrix as
				      * a square matrix with dimension @p{n}.
				      *
				      * In order to avoid the implicit
				      * conversion of integers and other types
				      * to a matrix, this constructor is
				      * declared @p{explicit}.
				      *
				      * By default, no memory is allocated.
				      */
    explicit FullMatrix (const unsigned int n = 0);
    
				     /**
				      * Constructor. Initialize the matrix as
				      * a rectangular matrix.
				      */
    FullMatrix (const unsigned int rows, const unsigned int cols);
    
				     /** 
				      * Copy constructor. This constructor does
				      * a deep copy of the matrix. Therefore,
				      * it poses an efficiency problem.
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
				      * Comparison operator. Be careful with
				      * this thing, it may eat up huge amounts
				      * of computing time! It is most commonly
				      * used for internal consistency checks
				      * of programs.
				      */
    bool operator == (const FullMatrix<number> &) const;

				     /**
				      * Assignment operator.
				      * Copy all elements of @p{src}
				      * into the matrix. The size is
				      * adjusted if needed.
				      *
				      * We can't use the other, templatized
				      * version since if we don't declare
				      * this one, the compiler will happily
				      * generate a predefined copy
				      * operator which is not what we want.
				      */
//    FullMatrix<number>& operator = (const FullMatrix<number>& src);

				     /**
				      * Assignment operator.
				      * Copy all elements of @p{src}
				      * into the matrix. The size is
				      * adjusted if needed.
				      */
//    template<typename number2>
//    FullMatrix<number>& operator = (const FullMatrix<number2>& src);


				     /**
				      * Fill rectangular block.
				      *
				      * The matrix @p{src} is copied
				      * into the target. The optional
				      * values @p{i} and @p{j} determine the
				      * upper left corner of the image
				      * of @p{src}.
				      *
				      * This function requires that
				      * @p{i+src.m()<=m()} and
				      * @p{j+src.n()<=n()}, that is, the
				      * image fits into the space of @p{this}.
				      */
    template<typename number2>
    void fill (const FullMatrix<number2> &src,
	       const unsigned int         i=0,
	       const unsigned int         j=0);
    

				     /**
				      * Make function of base class available.
				      */
    template<typename number2>
      void fill (const number2*);
    
				     /**
				      * Fill with permutation of another matrix.
				      *
				      * The matrix @p{src} is copied
				      * into the target. The two
				      * permutation @p{p_r} and
				      * @p{p_c} operate in a way, such
				      * that @{result(i,j) =
				      * src(p_r[i], p_c[j])}.
				      *
				      * The vectors may also be a
				      * selection from a larger set of
				      * integers, if the matrix
				      * @p{src} is bigger. It is also
				      * possible to duplicate rows or
				      * columns by this method.
				      */
    template<typename number2>
    void fill_permutation (const FullMatrix<number2> &src,
	       const std::vector<unsigned int>& p_rows,
	       const std::vector<unsigned int>& p_cols);

				     /**
				      * Set dimension to $m\times n$ and
				      * allocate memory if necessary. Forget
				      * the previous content of the matrix.
				      */
//    void reinit (const unsigned int m,
//		 const unsigned int n);
    
				     /**
				      * Set dimension to $n\times n$ and
				      * allocate memory if necessary. Forget
				      * the previous content of the matrix.
				      */
//    void reinit (const unsigned int n);
    
				     /**
				      * Set dimension to $m(B)\times n(B)$ and
				      * allocate memory if necessary. Forget
				      * the previous content of the matrix.
				      * However, this function does not copy
				      * the contents of @p{B}.
				      */
//    template<typename number2>
//    void reinit (const FullMatrix<number2> &B);
    
				     /**
				      * Number of rows of this matrix.
				      * To remember: this matrix is an
				      * $m \times n$-matrix.
				      */
    unsigned int m () const;
    
				     /**
				      * Number of columns of this matrix.
				      * To remember: this matrix is an
				      * $m \times n$-matrix.
				      */
    unsigned int n () const;

    				     /**
				      * Return whether the matrix contains only
				      * elements with value zero. This function
				      * is mainly for internal consistency
				      * checks and should seldomly be used when
				      * not in debug mode since it uses quite
				      * some time.
				      */
    bool all_zero () const;

    
				     /**
				      * Weighted addition. The matrix
				      * @p{s*B} is added to @p{this}.
				      *
				      * $A += sB$
				      */
    template<typename number2>
    void add (const number               s,
	      const FullMatrix<number2> &B);

				     /**
				      * Weighted addition of the
				      * transpose of @p{B} to @p{this}.
				      *
				      * $A += s B^T$
				      */
    template<typename number2>
    void Tadd (const number               s,
	       const FullMatrix<number2> &B);
    
				     /**
				      * Matrix-matrix-multiplication.
				      *
				      * The optional parameter
				      * @p{adding} determines, whether the
				      * result is stored in @p{C} or added
				      * to @p{C}.
				      *
				      * if (adding)
				      *  $C += A*B$
				      *
				      * if (!adding)
				      *  $C = A*B$
				      *
				      * Assumes that @p{A} and @p{B} have
				      * compatible sizes and that @p{C}
				      * already has the right size.
				      */
    template<typename number2>
    void mmult (FullMatrix<number2>       &C,
		const FullMatrix<number2> &B,
		const bool                 adding=false) const;
    
				     /**
				      * Matrix-matrix-multiplication using
				      * transpose of @p{this}.
				      *
				      * The optional parameter
				      * @p{adding} determines, whether the
				      * result is stored in @p{C} or added
				      * to @p{C}.
				      *
				      * if (adding)
				      *  $C += A^T*B$
				      *
				      * if (!adding)
				      *  $C = A^T*B$
				      *
				      * Assumes that @p{A} and @p{B} have
				      * compatible sizes and that @p{C}
				      * already has the right size.
				      */
    template<typename number2>
    void Tmmult (FullMatrix<number2>       &C,
		 const FullMatrix<number2> &B,
		 const bool                 adding=false) const;
    
				     /**
				      * Matrix-vector-multiplication.
				      *
				      * The optional parameter
				      * @p{adding} determines, whether the
				      * result is stored in @p{w} or added
				      * to @p{w}.
				      *
				      * if (adding)
				      *  $w += A*v$
				      *
				      * if (!adding)
				      *  $w = A*v$
				      */
    template<typename number2>
    void vmult (Vector<number2>       &w,
		const Vector<number2> &v,
		const bool             adding=false) const;
    
				     /**
				      * Transpose matrix-vector-multiplication.
				      * See @p{vmult} above.
				      */
    template<typename number2>
    void Tvmult (Vector<number2>       &w,
		 const Vector<number2> &v,
		 const bool             adding=false) const;

				     /**
				      * Return the square of the norm
				      * of the vector @p{v} with respect
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
				      * Obviously, the matrix needs to
				      * be square for this operation.
				      */
    template<typename number2>
    number2 matrix_norm_square (const Vector<number2> &v) const;

				     /**
				      * Build the matrix scalar product
				      * @p{u^T M v}. This function is mostly
				      * useful when building the cellwise
				      * scalar product of two functions in
				      * the finite element context.
				      */
    template<typename number2>
    number2 matrix_scalar_product (const Vector<number2> &u,
				   const Vector<number2> &v) const;

    				     /**
				      * Return the $l_1$-norm of the matrix, i.e.
				      * $|M|_1=max_{all columns j}\sum_{all 
				      * rows i} |M_ij|$,
				      * (max. sum of columns). This is the
				      * natural matrix norm that is compatible
				      * to the $l_1$-norm for vectors, i.e.
				      * $|Mv|_1\leq |M|_1 |v|_1$.
				      * (cf. Rannacher Numerik0)
				      */
    number l1_norm () const;

    				     /**
				      * Return the $l_\infty$-norm of the
				      * matrix, i.e.
				      * $|M|_\infty=\max_{all rows i}\sum_{all 
				      * columns j} |M_{ij}|$,
				      * (max. sum of rows).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the $l_\infty$-norm of vectors, i.e.
				      * $|Mv|_\infty \leq |M|_\infty |v|_\infty$.
				      * (cf. Rannacher Numerik0)
				      */
    number linfty_norm () const;
    
				     /**
				      * Compute the quadratic matrix norm.
				      * Return value is the root of the square
				      * sum of all matrix entries. Also called
				      * Frobenius norm.
				      * 
				      * This norm is compatible with the $l_2$
				      * vector norm. But it is not a natural
				      * matrix norm (cf Rannacher Numeric0),
				      * therefore it is not called $l_2$-norm.
				      */
    number norm2 () const;

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
    number relative_symmetry_norm2 () const;
    
				     /**
				      * A=Inverse(A). Inversion of
				      * this matrix by
				      * Gauss-Jordan algorithm with
				      * partial pivoting.  This
				      * process is well-behaved for
				      * positive definite matrices,
				      * but be aware of round-off
				      * errors in the indefinite case.
				      *
				      * The numerical effort to invert
				      * an @p{n x n} matrix is of the
				      * order @p{n**3}.
				      */
    void gauss_jordan ();

				     /**
                                      * Computes the determinant of a matrix.
                                      * This is only implemented for one two and
                                      * three dimensions, since for higher
                                      * dimensions the numerical work explodes.
                                      * Obviously, the matrix needs to be square
                                      * for this function.
                                      */
    double determinant () const;

				     /**
				      * Assign the inverse of the
				      * given matrix to
				      * @p{*this}. This function is
				      * hardcoded for square matrices
				      * of dimension one to four,
				      * since the amount of code
				      * needed grows quickly.  For
				      * larger matrices, the method
				      * @ref{gauss_jordan} is invoked
				      * implicitly.
				      */
    void invert (const FullMatrix<number> &M);

				     /**
				      * $A(i,1-n)+=s*A(j,1-n)$.
				      * Simple addition of rows of this
				      */
    void add_row (const unsigned int i,
		  const number       s,
		  const unsigned int j);

				     /**
				      * $A(i,1-n)+=s*A(j,1-n)+t*A(k,1-n)$.
				      * Multiple addition of rows of this.
				      */
    void add_row (const unsigned int i,
		  const number s, const unsigned int j,
		  const number t, const unsigned int k);

				     /**
				      * $A(1-n,i)+=s*A(1-n,j)$.
				      *  Simple addition of columns of this.
				      */
    void add_col (const unsigned int i,
		  const number       s,
		  const unsigned int j);

				     /**
				      * $A(1-n,i)+=s*A(1-n,j)+t*A(1-n,k)$.
				      *  Multiple addition of columns of this.
				      */
    void add_col (const unsigned int i,
		  const number s, const unsigned int j,
		  const number t, const unsigned int k);

				     /**
				      * Swap  A(i,1-n) <-> A(j,1-n).
				      * Swap rows i and j of this
				      */
    void swap_row (const unsigned int i, const unsigned int j);

				     /**
				      *  Swap  A(1-n,i) <-> A(1-n,j).
				      *  Swap columns i and j of this
				      */
    void swap_col (const unsigned int i, const unsigned int j);

				     /**
				      * $w=b-A*v$.
				      * Residual calculation , returns
				      * the $l_2$-norm $|w|$
				      */
    template<typename number2, typename number3>
    double residual (Vector<number2>      & w,
		     const Vector<number2>& v,
		     const Vector<number3>& b) const;

				     /**
				      * Forward elimination of lower triangle.
				      * Inverts the lower triangle of a
				      * quadratic matrix.
				      *
				      * If the matrix has more columns than rows,
				      * this function only operates on the left
				      * square submatrix. If there are more rows,
				      * the upper square part of the matrix
				      * is considered
				      */
    template<typename number2>
    void forward (Vector<number2>       &dst,
		  const Vector<number2> &src) const;

				     /**
				      * Backward elimination of upper triangle.
				      * @see forward
				      */
    template<typename number2>
    void backward (Vector<number2>       &dst,
		   const Vector<number2> &src) const;

				     /**
				      * QR-factorization of a matrix.
				      * The orthogonal transformation Q is
				      * applied to the vector y and this matrix.
				      *
				      * After execution of householder, the upper
				      * triangle contains the resulting matrix R,
				      * the lower the incomplete factorization
				      * matrices.
				      */
    template<typename number2>
    void householder (Vector<number2> &y);

				     /**
				      * Least-Squares-Approximation by
				      * QR-factorization.
				      */
    template<typename number2>
    double least_squares (Vector<number2> &dst,
			  Vector<number2> &src);

				     /**
				      *  A(i,i)+=B(i,1-n). Addition of complete
				      *  rows of B to diagonal-elements of this ; <p>
				      *  ( i = 1 ... m )
				      */
    template<typename number2>
    void add_diag (const number               s,
		   const FullMatrix<number2> &B);

				     /**
				      *  A(i,i)+=s  i=1-m.
				      * Add constant to diagonal elements of this
				      */
    void diagadd (const number s);

				     /**
				      * Output of the matrix in
				      * user-defined format.
				      */
    void print (std::ostream       &s,
		const unsigned int  width=5,
		const unsigned int  precision=2) const;

				     /**
				      * Print the matrix in the usual
				      * format, i.e. as a matrix and
				      * not as a list of nonzero
				      * elements. For better
				      * readability, elements not in
				      * the matrix are displayed as
				      * empty space, while matrix
				      * elements which are explicitely
				      * set to zero are displayed as
				      * such.
				      *
				      * The parameters allow for a
				      * flexible setting of the output
				      * format: @p{precision} and
				      * @p{scientific} are used to
				      * determine the number format,
				      * where @p{scientific} = @p{false}
				      * means fixed point notation.  A
				      * zero entry for @p{width} makes
				      * the function compute a width,
				      * but it may be changed to a
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
				      * This function
				      * may produce @em{large} amounts of
				      * output if applied to a large matrix!
				      */
    void print_formatted (std::ostream       &out,
			  const unsigned int  presicion=3,
			  const bool          scientific  = true,
			  const unsigned int  width       = 0,
			  const char         *zero_string = " ",
			  const double        denominator = 1.) const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyMatrix);
				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionMismatch,
		    int, int,
		    << "The two dimensions " << arg1 << " and " << arg2
		    << " do not match here.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotQuadratic);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotRegular);
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
    DeclException1 (ExcNotImplemented,
		    int,
		    << "This function is not implemented for the given"
		    << " matrix dimension " << arg1);

  private:
};


/*-------------------------Inline functions -------------------------------*/




template <typename number>
inline unsigned int
FullMatrix<number>::m() const
{
  return n_rows();
};



template <typename number>
inline unsigned int
FullMatrix<number>::n() const
{
  return n_cols();
};


template <typename number>
template <typename number2>
inline
void FullMatrix<number>::fill (const number2* src)
{
  vector2d<number>::fill(src);
}



/*----------------------------   fullmatrix.h     ---------------------------*/

#endif
/*----------------------------   fullmatrix.h     ---------------------------*/
