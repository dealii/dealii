//----------------------  sparse_decomposition.h  ---------------------------
//    Copyright (C) 1998, 1999, 2000, 2001, 2002
//    by the deal.II authors and Stephen "Cheffo" Kolaroff
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------  sparse_decomposition.h  ---------------------------
#ifndef __deal2__sparse_decomposition_h
#define __deal2__sparse_decomposition_h

#include <base/config.h>
#include <lac/sparse_matrix.h>

#include <cmath>

/**
 * Abstract base class for sparse LU decompositions of a sparse matrix
 * into another sparse matrix.
 *
 * The decomposition is stored as a sparse matrix, for
 * which the user has to give a sparsity pattern and which is why this
 * class is derived from the @p{SparseMatrix}. Since it is not a matrix in
 * the usual sense, the derivation is @p{protected} rather than @p{public}.
 *
 * @sect3{Fill-in}
 *
 * The sparse LU decompositions are frequently used with additional fill-in, i.e. the
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
 * sequence, for example (@p{lu_sparsity} is some sparsity pattern to be used
 * for the decomposition, which you have to create beforehand):
 * @begin{verbatim}
 *   SparseLUImplementation<double> lu (lu_sparsity);
 *   lu.decompose (global_matrix);
 *
 *   somesolver.solve (A, x, f,  lu);
 * @end{verbatim}
 *
 * @sect2{State management}
 *
 * In order to prevent users from applying decompositions before the
 * decomposition itself has been built, and to introduce some
 * optimization of common "sparse idioms", this class introduces a
 * simple state management.  A SparseLUdecomposition instance is
 * considered @p{not decomposed} if the decompose method has not yet
 * been invoked since the last time the underlying @ref{SparseMatrix}
 * had changed. The underlying sparse matrix is considered changed
 * when one of this class reinit methods, constructors or destructors
 * are invoked.  The @p{not decomposed} state is indicated by a false
 * value returned by @p{is_decomposed} method.  It is illegal to apply
 * this decomposition (@p{vmult} method) in not decomposed state; in
 * this case, the @p{vmult} method throws an @p{ExcInvalidState}
 * exception. This object turns into decomposed state immediately
 * after its @p{decompose} method is invoked. The @p{decomposed}
 * state is indicated by true value returned by @p{is_decomposed}
 * method. It is legal to apply this decomposition (@p{vmult} method) in
 * decomposed state.
 *
 *
 * @sect2{Particular implementations}
 *
 * It is enough to override the @p{decompose} and @p{vmult} methods to
 * implement particular LU decompositions, like the true LU, or the
 * Cholesky decomposition. Additionally, if that decomposition needs
 * fine tuned diagonal strengthening on a per row basis, it may override the
 * @p{get_strengthen_diagonal} method. You should invoke the non-abstract
 * base class method to employ the state management. Implementations
 * may choose more restrictive definition of what is legal or illegal
 * state; but they must conform to the @p{is_decomposed} method
 * specification above.
 *
 * If an exception is thrown by method other than @p{vmult}, this
 * object may be left in an inconsistent state.
 *
 * @author Stephen "Cheffo" Kolaroff, 2002, based on SparseILU implementation by Wolfgang Bangerth
 */
template <typename number>
class SparseLUDecomposition : protected SparseMatrix<number>
{
  public:

    				     /**
				      * Constructor; initializes the
				      * decomposition to be empty,
				      * without any structure, i.e.
				      * it is not usable at all. This
				      * constructor is therefore only
				      * useful for objects which are
				      * members of a class. All other
				      * matrices should be created at
				      * a point in the data flow where
				      * all necessary information is
				      * available.
				      *
				      * You have to initialize the
				      * matrix before usage with
				      * @p{reinit(SparsityPattern)}.
				      */

    SparseLUDecomposition ();

    				     /**
				      * Constructor. Takes the given
				      * matrix sparsity structure to
				      * represent the sparsity pattern
				      * of this decomposition.  You
				      * can change the sparsity
				      * pattern later on by calling
				      * the @p{reinit} function.
				      *
				      * You have to make sure that the
				      * lifetime of the sparsity
				      * structure is at least as long
				      * as that of this object or as
				      * long as @p{reinit} is not
				      * called with a new sparsity
				      * structure.
				      */
    SparseLUDecomposition (const SparsityPattern& sparsity);

                                     /**
                                      * Destruction.
                                      */
    virtual ~SparseLUDecomposition ();

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
				      * Perform the sparse LU
				      * factorization of the given
				      * matrix. After this method
				      * invocation, and before
				      * consecutive reinit invocation
				      * this object is in decomposed
				      * state.
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
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                    strengthen_diagonal=0.);

                                     /**
                                      * Determines if this object is
                                      * in synch with the underlying
                                      * @ref{SparsityPattern}.
                                      */ 
    virtual bool is_decomposed () const;	

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    virtual unsigned int memory_consumption () const;

                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcInvalidStrengthening,
		    double,
		    << "The strengthening parameter " << arg1
		    << " is not greater or equal than zero!");

                                     /**
                                      * Exception. Indicates violation
                                      * of a @p{state rule}.
                                      */
    DeclException0 (ExcInvalidState);

  protected:
                                     /**
                                      * Copies the passed SparseMatrix
                                      * onto this object. This
                                      * object's sparsity pattern
                                      * remains unchanged.
                                      */
    template<typename somenumber>
    void copy_from (const SparseMatrix<somenumber>& matrix);

                                     /**
                                      * Performs the strengthening
                                      * loop. For each row calculates
                                      * the sum of absolute values of
                                      * its elements, determines the
                                      * strengthening factor (through
                                      * @p{get_strengthen_diagonal})
                                      * sf and multiplies the diagonal
                                      * entry with @p{sf+1}.
                                      */
    virtual void strengthen_diagonal_impl ();

                                     /**
                                      * In the decomposition phase,
                                      * computes a strengthening
                                      * factor for the diagonal entry
                                      * in row @p{row} with sum of
                                      * absolute values of its
                                      * elements @p{rowsum}.<br> Note:
                                      * The default implementation in
                                      * @ref{SparseLUDecomposition}
                                      * returns
                                      * @p{strengthen_diagonal}'s
                                      * value.
                                      */
    virtual number get_strengthen_diagonal(const number rowsum, const unsigned int row) const;
    
                                     /**
                                      * State flag. If not in
                                      * @em{decomposed} state, it is
                                      * unlegal to apply the
                                      * decomposition.  This flag is
                                      * cleared when the underlaying
                                      * @ref{SparseMatrix}
                                      * @ref{SparsityPattern} is
                                      * changed, and set by
                                      * @p{decompose}.
                                      */
    bool decomposed;

                                     /**
                                      * The default strenghtening
                                      * value, returned by
                                      * @p{get_strengthen_diagonal}.
                                      */
    double  strengthen_diagonal;

                                     /**
                                      * For every row in the
                                      * underlying
                                      * @ref{SparsityPattern}, this
                                      * array contains a pointer
                                      * to the row's first
                                      * afterdiagonal entry. Becomes
                                      * available after invocation of
                                      * @p{decompose}.
                                      */
    std::vector<const unsigned int*> prebuilt_lower_bound;
    
  private:
                                     /**
                                      * Fills the
                                      * @ref{prebuilt_lower_bound}
                                      * array.
                                      */
    void prebuild_lower_bound ();
    
};



template <typename number>
inline number
SparseLUDecomposition<number>::
get_strengthen_diagonal(const number /*rowsum*/,
			const unsigned int /*row*/) const
{
  return strengthen_diagonal;
}



template <typename number>
inline bool
SparseLUDecomposition<number>::is_decomposed () const
{
  return decomposed;
}

#endif // __deal2__sparse_decomposition_h
