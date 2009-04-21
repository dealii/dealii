//---------------------------------------------------------------------------
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009 by the deal.II authors
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__sparse_ilu_h
#define __deal2__sparse_ilu_h


#include <base/config.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_decomposition.h>
#include <lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Incomplete LU decomposition of a sparse matrix into another sparse matrix.
 * A given matrix is decomposed into a incomplete LU factorization, where
 * by incomplete we mean that also a sparse decomposition is used and entries
 * in the decomposition that do not fit into the sparsity structure of this
 * object are discarded.
 *
 * The algorithm used by this class is essentially a copy of the
 * algorithm given in the book Y. Saad: "Iterative methods for sparse
 * linear systems", second edition, in section 10.3.2.
 *
 * 
 * <h3>Usage and state management</h3>
 *
 * Refer to SparseLUDecomposition documentation for suggested
 * usage and state management. This class is used in the @ref
 * step_22 "step-22" tutorial program.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 * 
 * @author Wolfgang Bangerth, 2008, 2009; unified interface: Ralf Hartmann
 */
template <typename number>
class SparseILU : public SparseLUDecomposition<number>
{
  public:
                                     /**
                                      * Constructor. Does nothing.
				      *
				      * Call the @p initialize
				      * function before using this
				      * object as preconditioner.
                                      */
    SparseILU ();

                                     /**
				      * @deprecated This method is
				      * deprecated, and left for
				      * backward compability. It will
				      * be removed in later versions.
                                      */
    SparseILU (const SparsityPattern &sparsity);

				     /**
				      * Make
				      * SparseLUDecomposition::AdditionalData
				      * accessible to this class as
				      * well.
				      */
    typedef
    typename SparseLUDecomposition<number>::AdditionalData
    AdditionalData;

				     /**
				      * Perform the incomplete LU
				      * factorization of the given
				      * matrix.
				      *
				      * This function needs to be
				      * called before an object of
				      * this class is used as
				      * preconditioner.
				      *
				      * For more details about
				      * possible parameters, see the
				      * class documentation of
				      * SparseLUDecomposition and the
				      * documentation of the
				      * @p SparseLUDecomposition::AdditionalData
				      * class.
				      *
				      * According to the
				      * @p parameters, this function
				      * creates a new SparsityPattern
				      * or keeps the previous sparsity
				      * or takes the sparsity given by
				      * the user to @p data. Then,
				      * this function performs the LU
				      * decomposition.
				      *
				      * After this function is called
				      * the preconditioner is ready to
				      * be used.
				      */
    template <typename somenumber>
    void initialize (const SparseMatrix<somenumber> &matrix,
		     const AdditionalData parameters = AdditionalData());

				     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be
				      * removed in later versions.
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                    strengthen_diagonal=0.);

				     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be
				      * removed in later versions.
				      */
    template <typename somenumber>
    void apply_decomposition (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src) const;

				     /**
				      * Apply the incomplete decomposition,
				      * i.e. do one forward-backward step
				      * $dst=(LU)^{-1}src$.
				      *
				      * The initialize() function
				      * needs to be called before.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>       &dst,
		const Vector<somenumber> &src) const;


				     /**
				      * Apply the transpose of the 
				      * incomplete decomposition,
				      * i.e. do one forward-backward step
				      * $dst=(LU)^{-T}src$.
				      *
				      * The initialize() function
				      * needs to be called before.
				      */
    template <typename somenumber>
    void Tvmult (Vector<somenumber>       &dst,
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
    DeclException1 (ExcInvalidStrengthening,
		    double,
		    << "The strengthening parameter " << arg1
		    << " is not greater or equal than zero!");
				     //@}
};

/*@}*/
//---------------------------------------------------------------------------

template <typename number>
template <typename somenumber>
inline
void
SparseILU<number>::apply_decomposition (Vector<somenumber>       &dst,
                                        const Vector<somenumber> &src) const
{
  vmult (dst, src);
}




DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_ilu_h
