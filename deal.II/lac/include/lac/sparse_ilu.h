//----------------------------  sparse_ilu.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_ilu.h  ---------------------------
#ifndef __deal2__sparse_ilu_h
#define __deal2__sparse_ilu_h


#include <lac/sparse_matrix.h>


/**
 * Incomplete LU decomposition of a sparse matrix into another sparse matrix.
 * A given matrix is decomposed into a incomplete LU factorization, where
 * by incomplete we mean that also a sparse decomposition is used and entries
 * in the decomposition that do not fit into the sparsity structure of this
 * object are discarded.
 *
 * The algorithm used by this class is as follows (indices run from @p{0}
 * to @p{N-1}):
 * \begin{verbatim}
 * copy original matrix into a[i,j]
 * 
 * for i=1..N-1
 *   a[i-1,i-1] = a[i-1,i-1]^{-1}
 *
 *   for k=0..i-1
 *     a[i,k] = a[i,k] * a[k,k]
 *
 *     for j=k+1..N-1
 *      if (a[i,j] exists & a[k,j] exists)
 *        a[i,j] -= a[i,k] * a[k,j]
 * \end{verbatim}
 * Using this algorithm, we store the decomposition as a sparse matrix, for
 * which the user has to give a sparsity pattern and which is why this
 * class is derived from the @p{SparseMatrix}. Since it is not a matrix in
 * the usual sense, the derivation is @p{protected} rather than @p{public}.
 *
 * Note that in the algorithm given, the lower left part of the matrix base
 * class is used to store the @p{L} part of the decomposition, while 
 * the upper right part is used to store @p{U}. The diagonal is used to
 * store the inverses of the diagonal elements of the decomposition; the
 * latter makes the application of the decomposition faster, since inversion
 * by the diagonal element has to be done only once, rather than at each
 * application (multiplication is much faster than division).
 *
 *
 * @sect3{Fill-in}
 *
 * The sparse ILU is frequently used with additional fill-in, i.e. the
 * sparsity structure of the decomposition is denser than that of the matrix
 * to be decomposed. The @p{decompose} function of this class allows this fill-in
 * as long as all entries present in the original matrix are present in the
 * decomposition also, i.e. the sparsity pattern of the decomposition is a
 * superset of the sparsity pattern in the original matrix.
 *
 * Such fill-in can be accomplished by various ways, one of which is a
 * copy-constructor of the @p{SparsityPattern} class which allows the addition
 * of side-diagonals to a given sparsity structure.
 *
 *
 * @sect3{Use as a preconditioner}
 *
 * If you want to use an object of this class as a preconditioner for another
 * matrix, you can do so by calling the solver function using the following
 * sequence, for example (@p{ilu_sparsity} is some sparsity pattern to be used
 * for the decomposition, which you have to create beforehand):
 * \begin{verbatim}
 *   SparseILU<double> ilu (ilu_sparsity);
 *   ilu.decompose (global_matrix);
 *
 *   somesolver.solve (A, x, f,
 *                     PreconditionUseMatrix<SparseILU<double>,Vector<double> >
 * 	                      (ilu,&SparseILU<double>::template apply_decomposition<double>));
 * \end{verbatim}
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
 * @author Wolfgang Bangerth, 1999, based on a similar implementation by Malte Braack
 */
template <typename number>
class SparseILU : protected SparseMatrix<number>
{
  public:
    				     /**
				      * Constructor; initializes the decomposition
				      * to be empty, without any structure, i.e.
				      * it is not usable at all. This
				      * constructor is therefore only useful
				      * for objects which are members of a
				      * class. All other matrices should be
				      * created at a point in the data flow
				      * where all necessary information is
				      * available.
				      *
				      * You have to initialize
				      * the matrix before usage with
				      * @p{reinit(SparsityPattern)}.
				      */
    SparseILU ();

    				     /**
				      * Constructor. Takes the given matrix
				      * sparsity structure to represent the
				      * sparsity pattern of this decomposition.
				      * You can change the sparsity pattern later
				      * on by calling the @p{reinit} function.
				      *
				      * You have to make sure that the lifetime
				      * of the sparsity structure is at least
				      * as long as that of this object or as
				      * long as @p{reinit} is not called with a
				      * new sparsity structure.
				      */
    SparseILU (const SparsityPattern &sparsity);

				     /**
				      * Reinitialize the object but keep to
				      * the sparsity pattern previously used.
				      * This may be necessary if you @p{reinit}'d
				      * the sparsity structure and want to
				      * update the size of the matrix.
				      *
				      * This function does nothing more than
				      * passing down to the sparse matrix
				      * object the call for the same function,
				      * which is necessary however, since that
				      * function is not publically visible
				      * any more.
				      */
    void reinit ();

				     /**
				      * Reinitialize the sparse matrix with the
				      * given sparsity pattern. The latter tells
				      * the matrix how many nonzero elements
				      * there need to be reserved.
				      *
				      * This function does nothing more than
				      * passing down to the sparse matrix
				      * object the call for the same function,
				      * which is necessary however, since that
				      * function is not publically visible
				      * any more.
				      */
    void reinit (const SparsityPattern &sparsity);

				     /**
				      * Perform the incomplete LU factorization
				      * of the given matrix.
				      *
				      * Note that the sparsity structures of
				      * the decomposition and the matrix passed
				      * to this function need not be equal,
				      * but that the pattern used by this
				      * matrix needs to contain all elements
				      * used by the matrix to be decomposed.
				      * Fill-in is thus allowed.
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                    strengthen_diagonal=0.);

				     /**
				      * Apply the incomplete decomposition,
				      * i.e. do one forward-backward step
				      * $dst=(LU)^{-1}src$.
				      */
    template <typename somenumber>
    void apply_decomposition (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src) const;

				     /**
				      * Same as @p{apply_decomposition}, format for LAC.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>       &dst,
		const Vector<somenumber> &src) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);
				     /**
				      * Exception
				      */
    DeclException2 (ExcSizeMismatch,
		    int, int,
		    << "The sizes " << arg1 << " and " << arg2
		    << " of the given objects do not match.");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidStrengthening,
		    double,
		    << "The strengthening parameter " << arg1
		    << " is not greater or equal than zero!");
};

template <typename number>
template <typename somenumber>
void
SparseILU<number>::vmult (Vector<somenumber>       &dst,
			  const Vector<somenumber> &src) const
{
  apply_decomposition(dst, src);
}


#endif
