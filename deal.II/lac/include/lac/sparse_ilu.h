//---------------------------------------------------------------------------
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
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
 * The algorithm used by this class is as follows (indices run from @p 0
 * to @p N-1):
 * @verbatim
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
 * @endverbatim
 *
 * 
 * @sect2{Usage and state management}
 *
 * Refer to SparseLUDecomposition documentation for suggested
 * usage and state management.
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
 * @author Wolfgang Bangerth, 1999, based on a similar implementation
 * by Malte Braack; unified interface: Ralf Hartmann
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
				      * object as preconditioner
				      * (@p vmult).
                                      */
    SparseILU ();

                                     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be
				      * removed in later versions.
                                      */
    SparseILU (const SparsityPattern &sparsity);

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
		     const AdditionalData parameters);

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
				      * The @p initialize function
				      * needs to be called beforehand.
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




#endif // __deal2__sparse_ilu_h
