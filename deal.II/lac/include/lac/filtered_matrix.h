//----------------------------  filtered_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  filtered_matrix.h  ---------------------------
#ifndef __deal2__filtered_matrix_h
#define __deal2__filtered_matrix_h



#include <base/smartpointer.h>
#include <base/thread_management.h>
#include <vector>
#include <algorithm>


template <typename number> class Vector;



/**
 * This class is a wrapper for linear systems of equations with simple
 * equality constraints fixing individual degrees of freedom to a
 * certain value such as when using Dirichlet boundary
 * values. Mathematically speaking, it is used to represent a system
 * of linear equations $Ax=b$ with the constraint that $B_D x = g_D$,
 * where $B_D$ is a rectangular matrix with exactly one $1$ in each
 * row, and these $1$s in those columns representing constrained
 * degrees of freedom (e.g. for Dirichlet boundary nodes, thus the
 * index $D$) and zeroes for all other diagonal entries, and $g_D$
 * having the requested nodal values for these constrained
 * nodes. Thus, the underdetermined equation $B_D x = g_D$ fixes only
 * the constrained nodes and does not impose any condition on the
 * others. We note that $B_D B_D^T = 1_D$, where $1_D$ is the identity
 * matrix with dimension as large as the number of constrained degrees
 * of freedom. Likewise, $B_D^T B_D$ is the diagonal matrix with
 * diagonal entries $0$ or $1$ that, when applied to a vector, leaves
 * all constrained nodes untouched and deletes all unconstrained ones.
 *
 * For solving such a system of equations, we first write down the
 * Lagrangian $L=1/2 x^T A x - x^T b + l^T B_D x$, where $l$
 * is a Lagrange multiplier for the constraints. The stationarity
 * condition then reads
 * \begin{verbatim}
 * [ A   B_D^T ] [x] = [b  ]
 * [ B_D 0     ] [l] = [g_D]
 * \begin{verbatim}
 *
 * The first equation then reads $B_D^T l = b-Ax$. On the other hand,
 * if we left-multiply the first equation by $B_D^T B_D$, we obtain
 * $B_D^T B_D A x + B_D^T l = B_D^T B_D b$ after equating $B_D B_D^T$
 * to the identity matrix. Inserting the previous equality, this
 * yields $(A - B_D^T B_D A) x = (1 - B_D^T B_D)b$. Since
 * $x=(1 - B_D^T B_D) x + B_D^T B_D x = (1 - B_D^T B_D) x + B_D^T g_D$,
 * we can restate the linear system:
 * $A_D x = (1 - B_D^T B_D)b - (1 - B_D^T B_D) A B^T g_D$, where
 * $A_D = (1 - B_D^T B_D) A (1 - B_D^T B_D)$ is the matrix where all
 * rows and columns corresponding to constrained nodes have been deleted.
 *
 * The last system of equation only defines the value of the
 * unconstrained nodes, while the constrained ones are determined by
 * the equation $B_D x = g_D$. We can combine these two linear systems
 * by using the zeroed out rows of $A_D$: if we set the diagonal to
 * $1$ and the corresponding zeroed out element of the right hand side
 * to that of $g_D$, then this fixes the constrained elements as
 * well. We can write this as follows:
 * $A_X x = (1 - B_D^T B_D)b - (1 - B_D^T B_D) A B^T g_D + B_D^T g_D$,
 * where $A_X = A_D + B_D^T B_D$. Note that the two parts of the
 * latter matrix operate on disjoint subspaces (the first on the
 * unconstrained nodes, the latter on the constrained ones).
 *
 * In iterative solvers, it is not actually necessary to compute $A_X$
 * explicitly, since only matrix-vector operations need to be
 * performed. This can be done in a three-step procedure that first
 * clears all elements in the incoming vector that belong to
 * constrained nodes, then performs the product with the matrix $A$,
 * then clears again. This class is a wrapper to this procedure, it
 * takes a pointer to a matrix with which to perform matrix-vector
 * products, and does the cleaning of constrained elements itself.
 * This class therefore implements an overloaded @p{vmult} function
 * that does the matrix-vector product, as well as @p{Tvmult} for
 * transpose matrix-vector multiplication and @p{residual} for
 * residual computation, and can thus be used as a matrix replacement
 * in lineaer solvers.
 *
 * It also has the ability to generate the modification of the right
 * hand side, through the @ref{apply_constraints} function.
 *
 *
 * @sect3{Connection to other classes}
 *
 * The function @p{MatrixTools::apply_boundary_values} does exactly
 * the same that this class does, except for the fact that that
 * function actually modifies the matrix. Due to this, it is only
 * possible to solve with a matrix onto which
 * @p{MatrixTools::apply_boundary_values} was applied for one right
 * hand side and one set of boundary values since the modification of
 * the right hand side depends on the original matrix.
 *
 * While this is fine (and the recommended way) in cases where only
 * one solution of the linear system is required, for example in
 * solving linear stationary systems, one would often like to have the
 * ability to solve multiply with the same matrix in nonlinear
 * problems (where one often does not want to update the Hessian
 * between Newton steps, despite having different right hand sides in
 * subsequent steps) or time dependent problems, without having to
 * re-assemble the matrix or copy it to temporary matrices with which
 * one then can work. For these cases, this class is meant.
 *
 *
 * @sect3{Usage}
 *
 * Usage is simple: create an object of this type, point it to a
 * matrix that shall be used for $A$ above (either through the
 * constructor, the copy constructor, or the
 * @ref{set_referenced_matrix} function), specify the list of boundary
 * values or other constraints (through the @ref{add_constraints}
 * function), and then for each required solution modify the right
 * hand side vector (through @ref{apply_constraints}) and use this
 * object as matrix object in a linear solver. As linear solvers
 * should only use @ref{vmult} and @ref{residual} functions of a
 * matrix class, this class should be as a good a matrix as any other
 * for that purpose.
 *
 * Furthermore, also the @ref{precondition_Jacobi} function is
 * provided (since the computation of diagonal elements of the
 * filtered matrix $A_X$ is simple), so you can use this as a
 * preconditioner. Some other function useful for matrices are also
 * available.
 *
 * A typical code snippet showing the above steps is as follows:
 * @begin{verbatim}
 *   ... // set up sparse matrix A and right hand side b somehow
 *
 *                     // initialize filtered matrix with
 *                     // matrix and boundary value constraints
 *   FilteredMatrix<SparseMatrix<double> > filtered_A (A);
 *   filtered_A.add_constraints (boundary_values);
 *
 *                     // set up a linear solver
 *   SolverControl control (1000, 1.e-10, false, false);
 *   PrimitiveVectorMemory<Vector<double> > mem;
 *   SolverCG<Vector<double> > solver (control, mem);
 *
 *                     // set up a preconditioner object
 *   PreconditionJacobi<FilteredMatrix<SparseMatrix<double> > > prec;
 *   prec.initialize (filtered_A, 1.2);
 *
 *                     // compute modification of right hand side
 *   filtered_A.apply_constraints (b, true);
 *
 *                     // solve for solution vector x
 *   solver.solve (filtered_A, x, b, prec);
 * @end{verbatim}
 *
 *
 * @sect3{Template arguments}
 *
 * This class takes as template arguments a matrix and a vector
 * class. The former must provide @p{vmult}, @p{Tvmult}, and
 * @p{residual} member function that operate on the vector type (the
 * second template argument). The latter template parameter must
 * provide access to indivual elements through @p{operator()},
 * assignment through @p{operator=}.
 *
 *
 * @sect3{Thread-safety}
 *
 * The functions that operate as a matrix and do not change the
 * internal state of this object are synchronised and thus
 * threadsafe. You need not serialize calls to @p{vmult} or
 * @p{residual} therefore. Because these functions require the use of
 * a temporary, they block mutual execution, however. It is necessary
 * to allocate this temporary vector in class space since otherwise we
 * would have to allocate such a vector each time one of the member
 * functions is called (which may be very often for slowly converging
 * linear systems), which would be a serious performance
 * bottleneck. If you don't want this serialization of operations, you
 * have to use several objects of this type.
 *
 * @author Wolfgang Bangerth 2001
 */
template <class Matrix, class Vector=::Vector<typename Matrix::value_type> >
class FilteredMatrix : public Subscriptor
{
  public:
				     /**
				      * Type of matrix entries. In
				      * analogy to the STL container
				      * classes.
				      */
    typedef typename Matrix::value_type value_type;

				     /**
				      * Typedef defining a type that
				      * represents a pair of degree of
				      * freedom index and the value it
				      * shall have.
				      */
    typedef typename std::pair<unsigned int,value_type> IndexValuePair;

				     /**
				      * Default constructor. You will
				      * have to set the matrix to be
				      * used later using the
				      * @{set_referenced_matrix}
				      * function.
				      */
    FilteredMatrix ();

				     /**
				      * Copy constructor. Use the
				      * matrix and the constraints set
				      * in the given object for the
				      * present one as well.
				      */
    FilteredMatrix (const FilteredMatrix &fm);
    
				     /**
				      * Constructor. Use the given
				      * matrix for future operations.
				      */
    FilteredMatrix (const Matrix &matrix);

				     /**
				      * Copy operator. Take over
				      * matrix and constraints from
				      * the other object.
				      */
    FilteredMatrix & operator = (const FilteredMatrix &fm);
    
				     /**
				      * Set the matrix to be used
				      * further on. You will probably
				      * also want to call the
				      * @ref{clear_constraints}
				      * function if constraits were
				      * previously added.
				      */
    void set_referenced_matrix (const Matrix &m);
    
				     /**
				      * Return a reference to the
				      * matrix that is used by this
				      * object.
				      */
    const Matrix & get_referenced_matrix () const;

				     /**
				      * Add a list of constraints to
				      * the ones already managed by
				      * this object. The actual data
				      * type of this list must be so
				      * that dereferenced iterators
				      * are pairs of indices and the
				      * corresponding values to be
				      * enforced on the respective
				      * solution vector's entry. Thus,
				      * the data type might be, for
				      * example, a @{std::list} or
				      * @p{std::vector} of
				      * @ref{IndexValuePair} objects,
				      * but also a
				      * @p{std::map<unsigned,value_type>}.
				      *
				      * It is an error if the argument
				      * contains an entry for a degree
				      * of freedom that has already
				      * been constrained
				      * previously. Furthermore, it is
				      * assumed that the list of
				      * constraints is sorted. If the
				      * input argument is a
				      * @p{std::map<unsigned,value_type>},
				      * this is automatically the
				      * case.
				      *
				      * Note that we do not check
				      * whether the input range is
				      * sorted, as this would be too
				      * expensive. You have to ensure
				      * this yourself.
				      */
    template <class ConstraintList>
    void add_constraints (const ConstraintList &new_constraints);

				     /**
				      * Delete the list of constraints
				      * presently in use.
				      */
    void clear_constraints ();

				     /**
				      * Apply the constraints to a
				      * right hand side vector. This
				      * needs to be done before
				      * starting to solve with the
				      * filtered matrix. If the matrix
				      * is symmetric, set the second
				      * parameter to @p{true} to use a
				      * faster algorithm.
				      */
    void apply_constraints (Vector     &v,
			    const bool  matrix_is_symmetric) const;

				     /**
				      * Return the dimension of the
				      * image space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int m () const;
    
				     /**
				      * Return the dimension of the
				      * range space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int n () const;

				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M*src$ with $M$
				      * being this matrix. (This
				      * matrix is the filtered one to
				      * which we store a reference.)
				      */
    void vmult (Vector       &dst,
		const Vector &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * @p{vmult} but takes the
				      * transposed matrix. (This
				      * matrix is the filtered one to
				      * which we store a reference.)
				      *
				      * Because we need to use a
				      * temporary variable and since
				      * we only allocate that each
				      * time the matrix changed, this
				      * function only works for square
				      * matrices.
				      */
    void Tvmult (Vector       &dst,
		 const Vector &src) const;
  
				     /**
				      * Return the square of the norm
				      * of the vector $v$ with respect
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
				      *
				      * Note that in many cases, you
				      * will not want to compute the
				      * norm with respect to the
				      * filtered matrix, but with
				      * respect to the original
				      * one. For example, if you want
				      * to compute the $L^2$ norm of a
				      * vector by forming the matrix
				      * norm with the mass matrix,
				      * then you want to use the
				      * original mass matrix, not the
				      * filtered one where you might
				      * have eliminated Dirichlet
				      * boundary values.
				      */
    value_type matrix_norm_square (const Vector &v) const;

				     /**
				      * Compute the residual of an
				      * equation @p{Mx=b}, where the
				      * residual is defined to be
				      * @p{r=b-Mx} with @p{x}
				      * typically being an approximate
				      * of the true solution of the
				      * equation. Write the residual
				      * into @p{dst}. The l2 norm of
				      * the residual vector is
				      * returned.
				      *
				      * Note that it is assumed that
				      * @{b} is a vector that has been
				      * treated by the
				      * @ref{modify_rhs} function,
				      * since we can then assume that
				      * the components of the residual
				      * which correspond to
				      * constrained degrees of freedom
				      * do not contribute to the
				      * residual at all.
				      */
    value_type residual (Vector       &dst,
			 const Vector &x,
			 const Vector &b) const;
    
				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the @p{src} vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * damping factor @p{omega}.
				      */
    void precondition_Jacobi (Vector           &dst,
			      const Vector     &src,
			      const value_type  omega = 1.) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object. Since we are
				      * not the owner of the matrix
				      * referenced, its memory
				      * consumption is not included.
				      */
    unsigned int memory_consumption () const;
    
  private:
				     /**
				      * Declare an abbreviation for an
				      * iterator into the array
				      * constraint pairs, since that
				      * data type is so often used and
				      * is rather awkward to write out
				      * each time.
				      */
    typedef typename std::vector<IndexValuePair>::const_iterator const_index_value_iterator;
    
				     /**
				      * Helper class used to sort
				      * pairs of indices and
				      * values. Only the index is
				      * considered as sort key.
				      */
    struct PairComparison 
    {
					 /**
					  * Function comparing the
					  * pairs @p{i1} and @p{i2}
					  * for their keys.
					  */
	bool operator () (const IndexValuePair &i1,
			  const IndexValuePair &i2) const;
    };
    
				     /**
				      * Pointer to the sparsity
				      * pattern used for this
				      * matrix. In order to guarantee
				      * that it is not deleted while
				      * still in use, we subscribe to
				      * it using the @p{SmartPointer}
				      * class.
				      */
    SmartPointer<const Matrix> matrix;

				     /**
				      * Sorted list of pairs denoting
				      * the index of the variable and
				      * the value to which it shall be
				      * fixed.
				      */
    typename std::vector<IndexValuePair> constraints;

				     /**
				      * Vector to be used as temporary
				      * storage. Since memory
				      * allocation is expensive, we do
				      * not want to allocate temporary
				      * vectors in each call to
				      * matrix-vector function, so we
				      * rather allocate it only once
				      * and then reuse it over and
				      * over again. Note that in a
				      * multithreaded environment, we
				      * have to synchronise access to
				      * this vector.
				      */
    mutable Vector tmp_vector;

				     /**
				      * Mutex used to synchronise use
				      * of the temporary vector.
				      */
    mutable Threads::ThreadMutex tmp_mutex;
    
				     /**
				      * Do the pre-filtering step,
				      * i.e. zero out those components
				      * that belong to constrained
				      * degrees of freedom.
				      */
    void pre_filter (Vector &v) const;

				     /**
				      * Do the postfiltering step,
				      * i.e. set constrained degrees
				      * of freedom to the value of the
				      * input vector, as the matrix
				      * contains only ones on the
				      * diagonal for these degrees of
				      * freedom.
				      */
    void post_filter (const Vector &in,
		      Vector       &out) const;

				     /**
				      * Based on the size of the
				      * matrix and type of the matrix
				      * and vector, allocate a
				      * temporary vector. This
				      * function has to be overloaded
				      * for the various template
				      * parameter choices.
				      */
    void allocate_tmp_vector ();

				     /**
				      * Determine all entries in the
				      * given column of the matrix
				      * except for the diagonal entry
				      * and return their index/value
				      * pairs. If the matrix is
				      * symmetric, use a faster
				      * algorithm.
				      *
				      * This function needs to be
				      * specialised for the different
				      * matrix types.
				      */
    void get_column_entries (const unsigned int           index,
			     typename std::vector<IndexValuePair> &column_entries,
			     const bool                   matrix_is_symmetric) const;
};


/*---------------------- Inline functions -----------------------------------*/


template <class Matrix, class Vector>
inline
bool
FilteredMatrix<Matrix,Vector>::PairComparison::
operator () (const IndexValuePair &i1,
	     const IndexValuePair &i2) const
{
  return (i1.first < i2.first);
};



template <class Matrix, class Vector>
template <class ConstraintList>
void
FilteredMatrix<Matrix,Vector>::
add_constraints (const ConstraintList &new_constraints)
{
				   // add new constraints to end
  const unsigned int old_size = constraints.size();
  constraints.reserve (old_size + new_constraints.size());
  constraints.insert (constraints.end(),
		      new_constraints.begin(),
		      new_constraints.end());
				   // then merge the two arrays to
				   // form one sorted one
  std::inplace_merge (constraints.begin(),
		      constraints.begin()+old_size,
		      constraints.end(),
		      PairComparison());
};



template <class Matrix, class Vector>
inline
const Matrix &
FilteredMatrix<Matrix,Vector>::get_referenced_matrix () const
{
  return *matrix;
};



template <class Matrix, class Vector>
inline
unsigned int FilteredMatrix<Matrix,Vector>::m () const
{
  return matrix->m();
};



template <class Matrix, class Vector>
inline
unsigned int FilteredMatrix<Matrix,Vector>::n () const
{
  return matrix->n();
};



/*----------------------------   filtered_matrix.h     ---------------------------*/

#endif
/*----------------------------   filtered_matrix.h     ---------------------------*/


