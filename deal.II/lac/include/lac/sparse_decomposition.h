//----------------------  sparse_decomposition.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004
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

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Abstract base class for sparse LU decompositions of a sparse matrix
 * into another sparse matrix.
 *
 * The decomposition is stored as a sparse matrix which is why this
 * class is derived from the SparseMatrix. Since it is not a matrix in
 * the usual sense, the derivation is <tt>protected</tt> rather than <tt>public</tt>.

 * @sect3{Fill-in}
 *
 * The sparse LU decompositions are frequently used with additional
 * fill-in, i.e. the sparsity structure of the decomposition is denser
 * than that of the matrix to be decomposed. The initialize()
 * function of this class allows this fill-in as long as all entries
 * present in the original matrix are present in the decomposition
 * also, i.e. the sparsity pattern of the decomposition is a superset
 * of the sparsity pattern in the original matrix.
 *
 * Such fill-in can be accomplished by various ways, one of which is a
 * copy-constructor of the SparsityPattern class which allows the addition
 * of side-diagonals to a given sparsity structure.
 *
 * @sect3{Unified use of preconditioners}
 *
 * An object of this class can be used in the same form as all
 * PreconditionBlock preconditioners:
 * @code
 * SparseLUDecomposition<double> lu;
 * lu.initialize(matrix, SparseLUDecomposition<double>::AdditionalData(...));
 *
 * somesolver.solve (A, x, f, lu);
 * @endcode
 *
 * Through the AdditionalData object it is possible to specify
 * additional parameters of the LU decomposition.
 *
 * 1/ The matrix diagonals can be strengthened by adding
 * <tt>strengthen_diagonal</tt> times the sum of the absolute row entries
 * of each row to the respective diagonal entries. By default no
 * strengthening is performed.
 *
 * 2/ By default, each initialize() function call creates its own
 * sparsity. For that, it copies the sparsity of <tt>matrix</tt> and adds a
 * specific number of extra off diagonal entries specified by
 * <tt>extra_off_diagonals</tt>.
 *
 * 3/ By setting <tt>use_previous_sparsity=true</tt> the sparsity is not
 * recreated but the sparsity of the previous initialize() call is
 * reused (recycled). This might be useful when several linear
 * problems on the same sparsity need to solved, as for example
 * several Newton iteration steps on the same triangulation. The
 * default is <tt>false</tt>.
 *
 * 4/ It is possible to give a user defined sparsity to
 * <tt>use_this_sparsity</tt>. Then, no sparsity is created but
 * <tt>*use_this_sparsity</tt> is used to store the decomposed matrix. For
 * restrictions on the sparsity see section `Fill-in' above).
 *
 *
 * @sect2{State management}
 *
 * The state management simply requires the initialize() function to
 * be called before the object is used as preconditioner.
 *
 * Obsolete:
 * In order to prevent users from applying decompositions before the
 * decomposition itself has been built, and to introduce some
 * optimization of common "sparse idioms", this class introduces a
 * simple state management.  A SparseLUdecomposition instance is
 * considered not decomposed if the decompose method has not yet
 * been invoked since the last time the underlying SparseMatrix
 * had changed. The underlying sparse matrix is considered changed
 * when one of this class reinit methods, constructors or destructors
 * are invoked.  The not-decomposed state is indicated by a false
 * value returned by is_decomposed() method.  It is illegal to apply
 * this decomposition (vmult() method) in not decomposed state; in
 * this case, the vmult() method throws an <tt>ExcInvalidState</tt>
 * exception. This object turns into decomposed state immediately
 * after its decompose() method is invoked. The decomposed
 * state is indicated by true value returned by is_decomposed()
 * method. It is legal to apply this decomposition (vmult() method) in
 * decomposed state.
 *
 * @sect2{Particular implementations}
 *
 * It is enough to override the initialize() and vmult() methods to
 * implement particular LU decompositions, like the true LU, or the
 * Cholesky decomposition. Additionally, if that decomposition needs
 * fine tuned diagonal strengthening on a per row basis, it may override the
 * get_strengthen_diagonal() method. You should invoke the non-abstract
 * base class method to employ the state management. Implementations
 * may choose more restrictive definition of what is legal or illegal
 * state; but they must conform to the is_decomposed() method
 * specification above.
 *
 * If an exception is thrown by method other than vmult(), this
 * object may be left in an inconsistent state.
 *
 * @author Stephen "Cheffo" Kolaroff, 2002, based on SparseILU implementation by Wolfgang Bangerth; unified interface: Ralf Hartmann, 2003
 */
template <typename number>
class SparseLUDecomposition : protected SparseMatrix<number>, 
                              public virtual Subscriptor
{
  public:

    				     /**
				      * Constructor. Does nothing.
				      *
				      * Call the initialize()
				      * function before using this
				      * object as preconditioner
				      * (vmult()).
				      */
    SparseLUDecomposition ();

    				     /**
				      * This method is deprecated, and
				      * left for backward
				      * compability. It will be removed
				      * in later versions.
				      */
    SparseLUDecomposition (const SparsityPattern& sparsity);

                                     /**
                                      * Destruction.
                                      */
    virtual ~SparseLUDecomposition ();

				     /**
				      * Deletes all member
				      * variables. Leaves the class in
				      * the state that it had directly
				      * after calling the constructor
				      */
    virtual void clear();
    
				     /**
				      * Parameters for
				      * SparseDecomposition.
				      */
    class AdditionalData
    {
      public:
					 /**
					  * Constructor. For the
					  * parameters' description,
					  * see below.
					  */
	AdditionalData (const double strengthen_diagonal=0,
			const unsigned int extra_off_diagonals=0,
			const bool use_previous_sparsity=false,
			const SparsityPattern *use_this_sparsity=0);

					 /**
					  * <tt>strengthen_diag</tt> times
					  * the sum of absolute row
					  * entries is added to the
					  * diagonal entries.
					  *
					  * Per default, this value is
					  * zero, i.e. the diagonal is
					  * not strengthened.
					  */
	double strengthen_diagonal;

					 /**
					  * By default, the
					  * <tt>initialize(matrix,
					  * data)</tt> function creates
					  * its own sparsity. This
					  * sparsity has the same
					  * SparsityPattern as
					  * <tt>matrix</tt> with some extra
					  * off diagonals the number
					  * of which is specified by
					  * <tt>extra_off_diagonals</tt>.
					  *
					  * The user can give a
					  * SparsityPattern to
					  * <tt>use_this_sparsity</tt>. Then
					  * this sparsity is used and
					  * the
					  * <tt>extra_off_diagonals</tt>
					  * argument is ignored.
					  */
	unsigned int extra_off_diagonals;

					 /**
					  * If this flag is true the
					  * initialize() function uses
					  * the same sparsity that was
					  * used during the previous
					  * initialize() call.
					  *
					  * This might be useful when
					  * several linear problems on
					  * the same sparsity need to
					  * solved, as for example
					  * several Newton iteration
					  * steps on the same
					  * triangulation.
					  */
	bool use_previous_sparsity;
	
					 /**
					  * When a
					  * SparsityPattern is
					  * given to this argument,
					  * the initialize()
					  * function calls
					  * <tt>reinit(*use_this_sparsity)</tt>
					  * causing this sparsity to
					  * be used.
					  *
					  * Note that the sparsity
					  * structures of
					  * <tt>*use_this_sparsity</tt> and
					  * the matrix passed to the
					  * initialize function need
					  * not be equal, but that the
					  * pattern used by this
					  * matrix needs to contain
					  * all elements used by the
					  * matrix to be decomposed.
					  * Fill-in is thus allowed.
					  */
	const SparsityPattern *use_this_sparsity;
    };

				     /**
				      * This function needs to be
				      * called before an object of
				      * this class is used as
				      * preconditioner.
				      *
				      * For more detail about possible
				      * parameters, see the class
				      * documentation and the
				      * documentation of the
				      * SparseLUDecomposition::AdditionalData
				      * class.
				      *
				      * According to the
				      * <tt>parameters</tt>, this function
				      * creates a new SparsityPattern
				      * or keeps the previous sparsity
				      * or takes the sparsity given by
				      * the user to <tt>data</tt>. Then,
				      * this function performs the LU
				      * decomposition.
				      *
				      * After this function is called
				      * the preconditioner is ready to
				      * be used (vmult).
				      */
    template <typename somenumber>
    void initialize (const SparseMatrix<somenumber> &matrix,
		     const AdditionalData parameters);

				     /**
				      * This method is deprecated,
				      * and left for backward
				      * compability. It will be removed
				      * in later versions.
				      */
    void reinit (const SparsityPattern &sparsity);

				     /**
				      * This method is deprecated,
				      * and left for backward
				      * compability. It will be removed
				      * in later versions.
				      */
    template <typename somenumber>
    void decompose (const SparseMatrix<somenumber> &matrix,
		    const double                    strengthen_diagonal=0.);

                                     /**
				      * This method is deprecated,
				      * and left for backward
				      * compability. It will be removed
				      * in later versions.
                                      */ 
    virtual bool is_decomposed () const;	

				     /**
				      * Return whether the object is
				      * empty. It calls the inherited
				      * SparseMatrix::empty() function.
				      */
    bool empty () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    virtual unsigned int memory_consumption () const;

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
                                      * get_strengthen_diagonal())
                                      * sf and multiplies the diagonal
                                      * entry with <tt>sf+1</tt>.
                                      */
    virtual void strengthen_diagonal_impl ();

                                     /**
                                      * In the decomposition phase,
                                      * computes a strengthening
                                      * factor for the diagonal entry
                                      * in row <tt>row</tt> with sum of
                                      * absolute values of its
                                      * elements <tt>rowsum</tt>.<br> Note:
                                      * The default implementation in
                                      * SparseLUDecomposition
                                      * returns
                                      * <tt>strengthen_diagonal</tt>'s
                                      * value.
                                      */
    virtual number get_strengthen_diagonal(const number rowsum, const unsigned int row) const;
    
                                     /**
                                      * State flag. If not in
                                      * decomposed state, it is
                                      * illegal to apply the
                                      * decomposition.  This flag is
                                      * cleared when the underlaying
                                      * SparseMatrix
                                      * SparsityPattern is
                                      * changed, and set by
                                      * decompose().
                                      */
    bool decomposed;

                                     /**
                                      * The default strenghtening
                                      * value, returned by
                                      * get_strengthen_diagonal().
                                      */
    double  strengthen_diagonal;

                                     /**
                                      * For every row in the
                                      * underlying
                                      * SparsityPattern, this
                                      * array contains a pointer
                                      * to the row's first
                                      * afterdiagonal entry. Becomes
                                      * available after invocation of
                                      * decompose().
                                      */
    std::vector<const unsigned int*> prebuilt_lower_bound;
    
  private:
                                     /**
                                      * Fills the
                                      * #prebuilt_lower_bound
                                      * array.
                                      */
    void prebuild_lower_bound ();

				     /**
				      * In general this pointer is
				      * zero except for the case that
				      * no SparsityPattern is
				      * given to this class. Then, a
				      * SparsityPattern is created
				      * and is passed down to the
				      * SparseMatrix base class.
				      *
				      * Nevertheless, the
				      * SparseLUDecomposition
				      * needs to keep ownership of
				      * this sparsity. It keeps this
				      * pointer to it enabling it to
				      * delete this sparsity at
				      * destruction time.
				      */
    SparsityPattern *own_sparsity;
};

/*@}*/
//----------------------------------------------------------------------//

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



template <typename number>
inline bool 
SparseLUDecomposition<number>::empty () const 
{
  return SparseMatrix<number>::empty(); 
}



//----------------------------------------------------------------------//


template <typename number>
SparseLUDecomposition<number>::AdditionalData::AdditionalData (
  const double strengthen_diag,
  const unsigned int extra_off_diag,
  const bool use_prev_sparsity,
  const SparsityPattern *use_this_spars):
		strengthen_diagonal(strengthen_diag),
		extra_off_diagonals(extra_off_diag),
		use_previous_sparsity(use_prev_sparsity),
		use_this_sparsity(use_this_spars)
{}


#endif // __deal2__sparse_decomposition_h
