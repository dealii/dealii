//----------------------------  dof_constraints.h  ---------------------------
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
//----------------------------  dof_constraints.h  ---------------------------
#ifndef __deal2__dof_constraints_h
#define __deal2__dof_constraints_h


#include <vector>
#include <utility>
#include <base/exceptions.h>
#include <base/subscriptor.h>

class SparsityPattern;
template <int rows, int columns> class BlockSparsityPattern;
template <typename number> class SparseMatrix;
template <typename number, int rows, int columns> class BlockSparseMatrix;
template <int n_blocks> class BlockIndices;

/**
 * This class represents the matrix denoting the distribution of the degrees
 * of freedom of hanging nodes.
 *
 * The matrix is organized in lines (rows), but only those lines are stored
 * where constraints are present. Lines where only one entry (identity) is
 * present are not stored if not explicitely inserted.
 *
 * Constraint matrices are used to handle hanging nodes and other constrained
 * degrees of freedom. When building the global system matrix and the right
 * hand sides, you normally build them without taking care of the constraints,
 * purely on a topological base, i.e. by a loop over cells. In order to do
 * actual calculations, you have to 'condense' these matrices: eliminate
 * constrained degrees of freedom and distribute the appropriate values to
 * the unconstrained dofs. This changes the sparsity pattern of the sparse
 * matrices used in finite element calculations und is thus a quite expensive
 * operation.
 *
 * Condensation is done in four steps: first the large matrix sparsity pattern
 * is created (e.g. using @ref{DoFHandler}@p{::create_sparsity_pattern}), then the
 * sparsity pattern of the condensed matrix is made out of the large sparsity
 * pattern and the constraints. After that the global matrix is assembled and
 * finally condensed. To do these steps, you have (at least) two possibilities:
 * @begin{itemize}
 * @item Use two different sparsity patterns and two different matrices: you
 *   may eliminate the lines and rows connected with a constraint and create
 *   a totally new sparsity pattern and a new system matrix. This has the
 *   advantage that the resulting system of equations is free from artifacts
 *   of the condensation process and is therefore faster in the solution process
 *   since no unnecessary multiplications occur (see below). However, there are
 *   two major drawbacks: keeping two matrices at the same time can be quite
 *   unacceptable in many cases, since these matrices may be several 10 or even
 *   100 MB large. Secondly, the condensation process is quite expensive, since
 *   @em{all} entries of the matrix have to be copied, not only those which are
 *   subject to constraints.
 *
 * @item Use only one sparsity pattern and one matrix: doing it this way, the
 *   condense functions add nonzero entries to the sparsity pattern of the
 *   large matrix (with constrained nodes in it) where the condensation process
 *   of the matrix will create additional nonzero elements. In the condensation
 *   process itself, lines and rows subject to constraints are distributed to
 *   the lines and rows of unconstrained nodes. The constrained lines remain in
 *   place, however, unlike in the first possibility described above. In order
 *   not to disturb the solution process, these lines and rows are filled with
 *   zeros and identity on the main diagonal; the appropriate value in the right
 *   hand sides is set to zero. This way, the constrained node will always get
 *   the value zero upon solution of the equation system and will not couple to
 *   other nodes any more.
 *
 *   This method has the advantage that only one matrix and sparsity pattern is
 *   needed thus using less memory. Additionally, the condensation process is
 *   less expensive, since not all but only constrained values in the matrix
 *   have to be copied. On the other hand, the solution process will take a bit
 *   longer, since matrix vector multiplications will incur multiplications
 *   with zeroes in the lines subject to constraints. Additionally, the vector
 *   size is larger than in the first possibility, resulting in more memory
 *   consumption for those iterative solution methods using a larger number of
 *   auxiliary vectors (e.g. methods using explicite orthogonalization
 *   procedures).
 * @end{itemize}
 *
 * Usually, the second way is chosen since memory consumption upon construction
 * of a second matrix rules out the first possibility.
 *
 * This class provides two sets of @p{condense} functions: those taking two
 * arguments refer to the first possibility above, those taking only one do
 * their job in-place and refer to the second possibility.
 *
 * Condensing vectors works exactly as described above for matrices.
 *
 * After solving the condensed system of equations, the solution vector has to
 * be redistributed. This is done by the two @p{distribute} function, one working
 * with two vectors, one working in-place. The operation of distribution undoes
 * the condensation process in some sense, but it should be noted that it is not
 * the inverse operation.
 *
 * @author Wolfgang Bangerth, 1998
 */
class ConstraintMatrix : public Subscriptor
{
  public:
				     /**
				      * Constructor
				      */
    ConstraintMatrix ();


				     /**
				      * Add a new line to the matrix.
				      */
    void add_line (const unsigned int line);

				     /**
				      * Add an entry to a given line. The list
				      * of lines is searched from the back to
				      * the front, so clever programming would
				      * add a new line (which is pushed to the
				      * back) and immediately afterwards fill
				      * the entries of that line. This way, no
				      * expensive searching is needed.
				      */
    void add_entry (const unsigned int line,
		    const unsigned int column,
		    const double value);

				     /**
				      * Close the filling of entries. Since the
				      * lines of a matrix of this type are
				      * usually filled in an arbitrary order and
				      * since we do not want to use associative
				      * constainers to store the lines, we need
				      * to sort the lines and within the lines
				      * the columns before usage of the matrix.
				      * This is done through this function.
				      *
				      * Also, zero entries are discarded, since
				      * they are not needed.
				      *
				      * After closing, no more entries are
				      * accepted.
				      */
    void close ();

				     /**
				      * Clear all entries of this matrix. Reset
				      * the flag determining whether new entries
				      * are accepted or not.
				      *
				      * This function may be called also on
				      * objects which are empty or already
				      * cleared.
				      */
    void clear ();

				     /**
				      * Return number of constraints stored in
				      * this matrix.
				      */
    unsigned int n_constraints () const;

				     /**
				      * Return whether the degree of
				      * freedom with number @p{index} is
				      * a constrained one.
				      *
				      * Note that if @p{close} was
				      * called before, then this
				      * function is significantly
				      * faster, since then the
				      * constrained degrees of freedom
				      * are sorted and we can do a
				      * binary search, while before
				      * @p{close} was called, we have to
				      * perform a linear search
				      * through all entries.
				      */
    bool is_constrained (const unsigned int index) const;

				     /**
				      * Return the maximum number of
				      * other dofs that one dof is
				      * constrained to. For example,
				      * in 2d a hanging node is
				      * constrained only to its two
				      * neighbors, so the returned
				      * value would be @p{2}. However,
				      * for higher order elements
				      * and/or higher dimensions, or
				      * other types of constraints,
				      * this number is no more
				      * obvious.
				      *
				      * The name indicates that within
				      * the system matrix, references
				      * to a constrained node are
				      * indirected to the nodes it is
				      * constrained to.
				      */
    unsigned int max_constraint_indirections () const;
    
				     /**
				      * Condense a given sparsity pattern. This
				      * function assumes the uncondensed
				      * matrix struct to be compressed and the
				      * one to be filled to be empty. The
				      * condensed structure is compressed
				      * afterwards.
				      *
				      * The constraint matrix object must be
				      * closed to call this function.
				      */
    void condense (const SparsityPattern &uncondensed,
		   SparsityPattern       &condensed) const;


				     /**
				      * This function does much the same as
				      * the above one, except that it condenses
				      * the matrix struct 'in-place'. It does
				      * not remove nonzero entries from the
				      * matrix but adds those needed for the
				      * process of distribution of the
				      * constrained degrees of freedom.
				      *
				      * Since this function adds new nonzero
				      * entries to the sparsity pattern, the
				      * argument must not be compressed. However
				      * the constraint matrix must be closed.
				      * The matrix struct is compressed at the
				      * end of the function.
				      */
    void condense (SparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square block sparsity
				      * patterns.
				      */
    template <int blocks>
    void condense (BlockSparsityPattern<blocks,blocks> &sparsity) const;
    
				     /**
				      * Condense a given matrix. The associated
				      * matrix struct should be condensed and
				      * compressed. It is the user's
				      * responsibility to guarantee that all
				      * entries in the @p{condensed} matrix be
				      * zero!
				      *
				      * The constraint matrix object must be
				      * closed to call this function.
				      */
    template<typename number>
    void condense (const SparseMatrix<number> &uncondensed,
		   SparseMatrix<number>       &condensed) const;

				     /**
				      * This function does much the same as
				      * the above one, except that it condenses
				      * the matrix 'in-place'. See the general
				      * documentation of this class for more
				      * detailed information.
				      */
    template<typename number>
    void condense (SparseMatrix<number> &matrix) const;

				     /**
				      * Same function as above, but
				      * condenses square block sparse
				      * matrices.
				      */
    template <typename number, int blocks>
    void condense (BlockSparseMatrix<number,blocks,blocks> &sparsity) const;
    
				     /**
				      * Condense the given vector @p{uncondensed}
				      * into @p{condensed}. It is the user's
				      * responsibility to guarantee that all
				      * entries of @p{condensed} be zero!
				      *
				      * The @p{VectorType} may be a
				      * @ref{Vector}@p{<float>},
				      * @ref{Vector}@p{<double>},
				      * @ref{BlockVector}@p{<...>}, or any
				      * other type having the same
				      * interface.
				      */
    template <class VectorType>
    void condense (const VectorType &uncondensed,
		   VectorType       &condensed) const;

				     /**
				      * Condense the given vector
				      * in-place. The @p{VectorType} may
				      * be a @ref{Vector}@p{<float>},
				      * @ref{Vector}@p{<double>},
				      * @ref{BlockVector}@p{<...>}, or any
				      * other type having the same
				      * interface.
				      */
    template <class VectorType>
    void condense (VectorType &vec) const;

				     /**
				      * Re-distribute the elements of
				      * the vector @p{condensed} to
				      * @p{uncondensed}. It is the
				      * user's responsibility to
				      * guarantee that all entries of
				      * @p{uncondensed} be zero!
				      *
				      * This function undoes the
				      * action of @p{condense} somehow,
				      * but it should be noted that it
				      * is not the inverse of
				      * @p{condense}.
				      *
				      * The @p{VectorType} may be a
				      * @ref{Vector}@p{<float>},
				      * @ref{Vector}@p{<double>},
				      * @ref{BlockVector}@p{<...>}, or any
				      * other type having the same
				      * interface.
				      */
    template <class VectorType>
    void distribute (const VectorType &condensed,
		     VectorType       &uncondensed) const;

				     /**
				      * Re-distribute the elements of
				      * the vector in-place. The
				      * @p{VectorType} may be a
				      * @ref{Vector}@p{<float>},
				      * @ref{Vector}@p{<double>},
				      * @ref{BlockVector}@p{<...>}, or any
				      * other type having the same
				      * interface.
				      */
    template <class VectorType>
    void distribute (VectorType &vec) const;
    
				     /**
				      * Delete hanging nodes in a
				      * vector.  Sets all hanging node
				      * values to zero. The
				      * @p{VectorType} may be a
				      * @ref{Vector}@p{<float>},
				      * @ref{Vector}@p{<double>},
				      * @ref{BlockVector}@p{<...>}, or any
				      * other type having the same
				      * interface.
				      */
    template <class VectorType>
    void set_zero (VectorType &vec) const;


				     /**
				      * Print the constraint lines. Mainly for
				      * debugging purposes.
				      *
				      * This function writes out all entries
				      * in the constraint matrix lines with
				      * their value in the form
				      * @p{row col : value}. Unconstrained lines
				      * containing only one identity entry are
				      * not stored in this object and are not
				      * printed.
				      */
    void print (ostream &) const;


				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixIsClosed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotClosed);
				     /**
				      * Exception
				      */
    DeclException1 (ExcLineInexistant,
		    unsigned int,
		    << "The specified line " << arg1
		    << " does not exist.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);
				     /**
				      * Exception
				      */
    DeclException0 (ExcWrongDimension);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /**
				      * Exception
				      */
    DeclException4 (ExcEntryAlreadyExists,
		    int, int, double, double,
		    << "The entry for the indices " << arg1 << " and "
		    << arg2 << " already exists, but the values "
		    << arg3 << " (old) and " << arg4 << " (new) differ.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcDoFConstrainedToConstrainedDoF,
		    int, int,
		    << "You tried to constrain DoF " << arg1
		    << " to DoF " << arg2
		    << ", but that one is also constrained. This is not allowed!");
    
  private:

				     /**
				      * This class represents one line of a
				      * constraint matrix.
				      */
    struct ConstraintLine
    {
					 /**
					  * Number of this line. Since only very
					  * few lines are stored, we can not
					  * assume a specific order and have
					  * to store the line number explicitely.
					  */
	unsigned int line;

					 /**
					  * Row numbers and values of the entries
					  * in this line.
					  *
					  * For the reason why we use a vector
					  * instead of a map and the consequences
					  * thereof, the same applies as what is
					  * said for @ref{ConstraintMatrix}@p{::lines}.
					  */
	vector<pair<unsigned int,double> > entries;

					 /**
					  * This operator is a bit
					  * weird and unintuitive: it
					  * compares the line numbers
					  * of two lines. We need this
					  * to sort the lines; in fact
					  * we could do this using a
					  * comparison predicate.
					  * However, this way, it is
					  * easier, albeit unintuitive
					  * since two lines really
					  * have no god-given order
					  * relation.
					  */
	bool operator < (const ConstraintLine &) const;
    };

				     /**
				      * Store the lines of the matrix.
				      * Entries are usually
				      * appended in an arbitrary order and
				      * insertion into a vector is done best
				      * at the end, so the order is
				      * unspecified after all entries are
				      * inserted. Sorting of the entries takes
				      * place when calling the @p{close()} function.
				      *
				      * We could, instead of using a vector, use
				      * an associative array, like a map to
				      * store the lines. This, however, would
				      * mean a much more fractioned heap since it
				      * allocates many small objects, ans would
				      * additionally make usage of this matrix
				      * much slower.
				      */
    vector<ConstraintLine> lines;
	
				     /**
				      * Store whether the arrays are sorted.
				      * If so, no new entries can be added.
				      */
    bool sorted;
};


#endif
