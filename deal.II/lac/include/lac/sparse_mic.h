//----------------------------  sparse_mic.h  ---------------------------
//    Copyright (C) 1998, 1999, 2000, 2001, 2002
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_mic.h  ---------------------------
#ifndef __deal2__sparse_mic_h
#define __deal2__sparse_mic_h

#include <lac/sparse_matrix.h>
#include <lac/sparse_decomposition.h>


/**
 * Modified incomplete Cholesky (MIC(0)) preconditioner.  This class
 * conforms to the state and usage specification in
 * @ref{SparseLUDecomposition}.
 *
 * 
 * @sect2{The decomposition}
 * 
 * Let a sparse matrix A is in the form A = - L - U + D, where -L and
 * -U are strictly lower and upper triangular matrices. The MIC(0)
 * decomposition of the matrix A is defined by B = (X-L)X^(-1)(X-U),
 * where X is a diagonal matrix, defined by the condition rowsum(A) =
 * rowsum(B).
 * 
 * @author Stephen "Cheffo" Kolaroff, 2002.
 */
template <typename number>
class SparseMIC : public SparseLUDecomposition<number>
{
  public:
                                     /**
                                      * Constructor. Does nothing, so
                                      * you have to call @p{reinit}
                                      * sometimes afterwards.
                                      */
    SparseMIC ();

                                     /**
                                      * Constructor. Initialize the
                                      * sparsity pattern of this
                                      * object with the given
                                      * argument.
                                      */
    SparseMIC (const SparsityPattern &sparsity);
    
				     /**
				      * Make the @p{AdditionalData}
				      * type in the base class
				      * accessible to this class as
				      * well.
				      */
    typedef
    typename SparseLUDecomposition<number>::AdditionalData
    AdditionalData;

				     /**
				      * Reinitialize the object but
				      * keep to the sparsity pattern
				      * previously used.  This may be
				      * necessary if you @p{reinit}'d
				      * the sparsity structure and
				      * want to update the size of the
				      * matrix.
				      *
				      * After this method is invoked,
				      * this object is out of synch
				      * (not decomposed state).
				      *
				      * This function only releases
				      * some memory and calls the
				      * respective function of the
				      * base class.
				      */
    void reinit ();

				     /**
				      * Reinitialize the sparse matrix
				      * with the given sparsity
				      * pattern. The latter tells the
				      * matrix how many nonzero
				      * elements there need to be
				      * reserved.
				      *
				      *
				      * This function only releases
				      * some memory and calls the
				      * respective function of the
				      * base class.
				      */
    void reinit (const SparsityPattern &sparsity);

				     /**
				      * Same as @p{decompose}.
				      */
    template <typename somenumber>
    void initialize (const SparseMatrix<somenumber> &matrix,
		     const AdditionalData parameters);

    				     /**
				      * Perform the incomplete LU
				      * factorization of the given
				      * matrix.
				      *
				      * Note that the sparsity
				      * structures of the
				      * decomposition and the matrix
				      * passed to this function need
				      * not be equal, but that the
				      * pattern used by this matrix
				      * needs to contain all elements
				      * used by the matrix to be
				      * decomposed.  Fill-in is thus
				      * allowed.
				      *
				      * If @p{strengthen_diagonal}
				      * parameter is greater than
				      * zero, this method invokes the
				      * @p{strengthen_diagonal_impl}
				      * function of the base class.
				      *
				      * Refer to
				      * @ref{SparseLUDecomposition}
				      * documentation for state
				      * management.
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                   strengthen_diagonal=0.);

				     /**
				      * Apply the incomplete decomposition,
				      * i.e. do one forward-backward step
				      * $dst=(LU)^{-1}src$.
				      *
				      * Refer to
				      * @ref{SparseLUDecomposition}
				      * documentation for state
				      * management.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>       &dst,
                const Vector<somenumber> &src) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcMatrixNotSquare);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcStrengthenDiagonalTooSmall);
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
                                     /**
                                      * Exception
                                      */
    DeclException2(ExcDecompositionNotStable, int, double,
		   << "The diagonal element (" <<arg1<<","<<arg1<<") is "
		   << arg2 <<", but must be positive");
    
  private:
                                     /**
                                      * Values of the computed
                                      * diagonal.
                                      */
    std::vector<number> diag;
    
                                     /**
                                      * Inverses of the the diagonal:
                                      * precomputed for faster vmult.
                                      */
    std::vector<number> inv_diag;

                                     /**
                                      * Values of the computed "inner
                                      * sums", i.e. per-row sums of
                                      * the elements laying on the
                                      * right side of the diagonal.
                                      */
    std::vector<number> inner_sums;
    
                                     /**
                                      * Compute the row-th "inner
                                      * sum".
                                      */
    number get_rowsum (const unsigned int row) const;
};



template <typename number>
template <typename somenumber>
inline
void SparseMIC<number>::initialize (const SparseMatrix<somenumber> &matrix,
				    const AdditionalData data)
{
  decompose(matrix, data.strengthen_diagonal);
}



#endif  // __deal2__
