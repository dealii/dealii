//---------------------------------------------------------------------------
//    Copyright (C) 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__sparse_mic_h
#define __deal2__sparse_mic_h

#include <lac/sparse_matrix.h>
#include <lac/sparse_decomposition.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Modified incomplete Cholesky (MIC(0)) preconditioner.  This class
 * conforms to the state and usage specification in
 * SparseLUDecomposition.
 *
 * 
 * <h3>The decomposition</h3>
 * 
 * Let a sparse matrix A is in the form A = - L - U + D, where -L and
 * -U are strictly lower and upper triangular matrices. The MIC(0)
 * decomposition of the matrix A is defined by B = (X-L)X^(-1)(X-U),
 * where X is a diagonal matrix, defined by the condition rowsum(A) =
 * rowsum(B).
 * 
 * @author Stephen "Cheffo" Kolaroff, 2002, unified interface: Ralf
 * Hartmann 2003.
 */
template <typename number>
class SparseMIC : public SparseLUDecomposition<number>
{
  public:
                                     /**
                                      * Constructor. Does nothing, so
                                      * you have to call @p reinit
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
                                      * Destruction.
                                      */
    virtual ~SparseMIC();

				     /**
				      * Deletes all member
				      * variables. Leaves the class in
				      * the state that it had directly
				      * after calling the constructor
				      */
    virtual void clear();

				     /**
				      * Make the @p AdditionalData
				      * type in the base class
				      * accessible to this class as
				      * well.
				      */
    typedef
    typename SparseLUDecomposition<number>::AdditionalData
    AdditionalData;

				     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be
				      * removed in later versions.
				      */
    void reinit (const SparsityPattern &sparsity);

				     /**
				      * Same as @p decompose.
				      */
    template <typename somenumber>
    void initialize (const SparseMatrix<somenumber> &matrix,
		     const AdditionalData parameters);

    				     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be
				      * removed in later versions.
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                   strengthen_diagonal=0.);

				     /**
				      * Apply the incomplete decomposition,
				      * i.e. do one forward-backward step
				      * $dst=(LU)^{-1}src$.
				      *
				      * Call @p initialize before
				      * calling this function.
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

    				     /** @addtogroup Exceptions
				      * @{ */

                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcStrengthenDiagonalTooSmall);
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

				     //@}
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

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif  // __deal2__
