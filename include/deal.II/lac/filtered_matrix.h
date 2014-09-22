// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__filtered_matrix_h
#define __deal2__filtered_matrix_h



#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector_memory.h>
#include <vector>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;
template <class VECTOR> class FilteredMatrixBlock;


/*! @addtogroup Matrix2
 *@{
 */


/**
 * This class is a wrapper for linear systems of equations with simple
 * equality constraints fixing individual degrees of freedom to a
 * certain value such as when using Dirichlet boundary
 * values.
 *
 * In order to accomplish this, the vmult(), Tvmult(), vmult_add() and
 * Tvmult_add functions modify the same function of the original
 * matrix such as if all constrained entries of the source vector were
 * zero. Additionally, all constrained entries of the destination
 * vector are set to zero.
 *
 * <h3>Usage</h3>
 *
 * Usage is simple: create an object of this type, point it to a
 * matrix that shall be used for $A$ above (either through the
 * constructor, the copy constructor, or the
 * set_referenced_matrix() function), specify the list of boundary
 * values or other constraints (through the add_constraints()
 * function), and then for each required solution modify the right
 * hand side vector (through apply_constraints()) and use this
 * object as matrix object in a linear solver. As linear solvers
 * should only use vmult() and residual() functions of a
 * matrix class, this class should be as good a matrix as any other
 * for that purpose.
 *
 * Furthermore, also the precondition_Jacobi() function is
 * provided (since the computation of diagonal elements of the
 * filtered matrix $A_X$ is simple), so you can use this as a
 * preconditioner. Some other functions useful for matrices are also
 * available.
 *
 * A typical code snippet showing the above steps is as follows:
 * @code
 *   ... // set up sparse matrix A and right hand side b somehow
 *
 *                     // initialize filtered matrix with
 *                     // matrix and boundary value constraints
 *   FilteredMatrix<Vector<double> > filtered_A (A);
 *   filtered_A.add_constraints (boundary_values);
 *
 *                     // set up a linear solver
 *   SolverControl control (1000, 1.e-10, false, false);
 *   GrowingVectorMemory<Vector<double> > mem;
 *   SolverCG<Vector<double> > solver (control, mem);
 *
 *                     // set up a preconditioner object
 *   PreconditionJacobi<SparseMatrix<double> > prec;
 *   prec.initialize (A, 1.2);
 *   FilteredMatrix<Vector<double> > filtered_prec (prec);
 *   filtered_prec.add_constraints (boundary_values);
 *
 *                     // compute modification of right hand side
 *   filtered_A.apply_constraints (b, true);
 *
 *                     // solve for solution vector x
 *   solver.solve (filtered_A, x, b, filtered_prec);
 * @endcode
 *
 *
 * <h3>Connection to other classes</h3>
 *
 * The function MatrixTools::apply_boundary_values() does exactly
 * the same that this class does, except for the fact that that
 * function actually modifies the matrix. Consequently, it is only
 * possible to solve with a matrix to which
 * MatrixTools::apply_boundary_values() was applied for one right
 * hand side and one set of boundary values since the modification of
 * the right hand side depends on the original matrix.
 *
 * While this is a feasible method in cases where only
 * one solution of the linear system is required, for example in
 * solving linear stationary systems, one would often like to have the
 * ability to solve multiple times with the same matrix in nonlinear
 * problems (where one often does not want to update the Hessian
 * between Newton steps, despite having different right hand sides in
 * subsequent steps) or time dependent problems, without having to
 * re-assemble the matrix or copy it to temporary matrices with which
 * one then can work. For these cases, this class is meant.
 *
 *
 * <h3>Some background</h3>
 * Mathematically speaking, it is used to represent a system
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
 * @code
 * [ A   B_D^T ] [x] = [b  ]
 * [ B_D 0     ] [l] = [g_D]
 * @endcode
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
 * This class therefore implements an overloaded @p vmult function
 * that does the matrix-vector product, as well as @p Tvmult for
 * transpose matrix-vector multiplication and @p residual for
 * residual computation, and can thus be used as a matrix replacement
 * in linear solvers.
 *
 * It also has the ability to generate the modification of the right
 * hand side, through the apply_constraints() function.
 *
 *
 *
 * <h3>Template arguments</h3>
 *
 * This class takes as template arguments a matrix and a vector
 * class. The former must provide @p vmult, @p vmult_add,  @p Tvmult, and
 * @p residual member function that operate on the vector type (the
 * second template argument). The latter template parameter must
 * provide access to indivual elements through <tt>operator()</tt>,
 * assignment through <tt>operator=</tt>.
 *
 *
 * <h3>Thread-safety</h3>
 *
 * The functions that operate as a matrix and do not change the internal state
 * of this object are synchronised and thus threadsafe. Consequently, you do
 * not need to serialize calls to @p vmult or @p residual .
 *
 * @author Wolfgang Bangerth 2001, Luca Heltai 2006, Guido Kanschat 2007, 2008
 */
template <class VECTOR>
class FilteredMatrix : public Subscriptor
{
public:
  class const_iterator;

  /**
   * Declare the type of container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Accessor class for iterators
   */
  class Accessor
  {
    /**
     * Constructor. Since we use
     * accessors only for read
     * access, a const matrix
     * pointer is sufficient.
     */
    Accessor (const FilteredMatrix<VECTOR> *matrix,
              const size_type               index);

  public:
    /**
     * Row number of the element
     * represented by this
     * object.
     */
    size_type row() const;

    /**
     * Column number of the
     * element represented by
     * this object.
     */
    size_type column() const;

    /**
     * Value of the right hand
     * side for this row.
     */
    double value() const;

  private:
    /**
     * Advance to next entry
     */
    void advance ();

    /**
     * The matrix accessed.
     */
    const FilteredMatrix<VECTOR> *matrix;

    /**
     * Current row number.
     */
    size_type index;
    /*
     * Make enclosing class a
     * friend.
     */
    friend class const_iterator;
  };

  /**
   * STL conforming iterator.
   */
  class const_iterator
  {
  public:
    /**
     * Constructor.
     */
    const_iterator(const FilteredMatrix<VECTOR> *matrix,
                   const size_type index);

    /**
     * Prefix increment.
     */
    const_iterator &operator++ ();

    /**
     * Postfix increment.
     */
    const_iterator &operator++ (int);

    /**
     * Dereferencing operator.
     */
    const Accessor &operator* () const;

    /**
     * Dereferencing operator.
     */
    const Accessor *operator-> () const;

    /**
     * Comparison. True, if
     * both iterators point to
     * the same matrix
     * position.
     */
    bool operator == (const const_iterator &) const;
    /**
     * Inverse of <tt>==</tt>.
     */
    bool operator != (const const_iterator &) const;

    /**
     * Comparison operator. Result is
     * true if either the first row
     * number is smaller or if the row
     * numbers are equal and the first
     * index is smaller.
     */
    bool operator < (const const_iterator &) const;

    /**
     * Comparison operator. Compares just
     * the other way around than the
     * operator above.
     */
    bool operator > (const const_iterator &) const;

  private:
    /**
     * Store an object of the
     * accessor class.
     */
    Accessor accessor;
  };

  /**
   * Typedef defining a type that
   * represents a pair of degree of
   * freedom index and the value it
   * shall have.
   */
  typedef std::pair<size_type, double> IndexValuePair;

  /**
   * @name Constructors and initialization
   */
//@{
  /**
   * Default constructor. You will
   * have to set the matrix to be
   * used later using
   * initialize().
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
   *
   * @arg @p m: The matrix being used in multiplications.
   *
   * @arg @p
   * expect_constrained_source: See
   * documentation of
   * #expect_constrained_source.
   */
  template <class MATRIX>
  FilteredMatrix (const MATRIX &matrix,
                  bool expect_constrained_source = false);

  /**
   * Copy operator. Take over
   * matrix and constraints from
   * the other object.
   */
  FilteredMatrix &operator = (const FilteredMatrix &fm);

  /**
   * Set the matrix to be used
   * further on. You will probably
   * also want to call the
   * clear_constraints()
   * function if constraits were
   * previously added.
   *
   * @arg @p m: The matrix being used in multiplications.
   *
   * @arg @p
   * expect_constrained_source: See
   * documentation of
   * #expect_constrained_source.
   */
  template <class MATRIX>
  void initialize (const MATRIX &m,
                   bool expect_constrained_source = false);

  /**
   * Delete all constraints and the
   * matrix pointer.
   */
  void clear ();
//@}
  /**
   * @name Managing constraints
   */
//@{
  /**
   * Add the constraint that the
   * value with index <tt>i</tt>
   * should have the value
   * <tt>v</tt>.
   */
  void add_constraint (const size_type i, const double v);

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
   * example, a @p std::list or
   * @p std::vector of
   * IndexValuePair objects,
   * but also a
   * <tt>std::map<unsigned, double></tt>.
   *
   * The second component of these
   * pairs will only be used in
   * apply_constraints(). The first
   * is used to set values to zero
   * in matrix vector
   * multiplications.
   *
   * It is an error if the argument
   * contains an entry for a degree
   * of freedom that has already
   * been constrained
   * previously.
   */
  template <class ConstraintList>
  void add_constraints (const ConstraintList &new_constraints);

  /**
   * Delete the list of constraints
   * presently in use.
   */
  void clear_constraints ();
//@}
  /**
   * Vector operations
   */
//@{
  /**
   * Apply the constraints to a
   * right hand side vector. This
   * needs to be done before
   * starting to solve with the
   * filtered matrix. If the matrix
   * is symmetric (i.e. the matrix
   * itself, not only its sparsity
   * pattern), set the second
   * parameter to @p true to use a
   * faster algorithm.
   */
  void apply_constraints (VECTOR     &v,
                          const bool  matrix_is_symmetric) const;

  /**
   * Matrix-vector multiplication:
   * this operation performs
   * pre_filter(), multiplication
   * with the stored matrix and
   * post_filter() in that order.
   */
  void vmult (VECTOR       &dst,
              const VECTOR &src) const;

  /**
   * Matrix-vector multiplication:
   * this operation performs
   * pre_filter(), transposed
   * multiplication with the stored
   * matrix and post_filter() in
   * that order.
   */
  void Tvmult (VECTOR       &dst,
               const VECTOR &src) const;

  /**
   * Adding matrix-vector multiplication.
   *
   * @note The result vector of
   * this multiplication will have
   * the constraint entries set to
   * zero, independent of the
   * previous value of
   * <tt>dst</tt>. We excpect that
   * in most cases this is the
   * required behavior.
   */
  void vmult_add (VECTOR       &dst,
                  const VECTOR &src) const;

  /**
   * Adding transpose matrix-vector multiplication:
   *
   * @note The result vector of
   * this multiplication will have
   * the constraint entries set to
   * zero, independent of the
   * previous value of
   * <tt>dst</tt>. We excpect that
   * in most cases this is the
   * required behavior.
   */
  void Tvmult_add (VECTOR       &dst,
                   const VECTOR &src) const;
//@}

  /**
   * @name Iterators
   */
//@{
  /**
   * Iterator to the first
   * constraint.
   */
  const_iterator begin () const;
  /**
   * Final iterator.
   */
  const_iterator end () const;
//@}

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object. Since we are
   * not the owner of the matrix
   * referenced, its memory
   * consumption is not included.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Determine, whether
   * multiplications can expect
   * that the source vector has all
   * constrained entries set to
   * zero.
   *
   * If so, the auxiliary vector
   * can be avoided and memory as
   * well as time can be saved.
   *
   * We expect this for instance in
   * Newton's method, where the
   * residual already should be
   * zero on constrained
   * nodes. This is, because there
   * is no testfunction in these
   * nodes.
   */
  bool expect_constrained_source;

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
     * pairs @p i1 and @p i2
     * for their keys.
     */
    bool operator () (const IndexValuePair &i1,
                      const IndexValuePair &i2) const;
  };

  /**
   * Pointer to the sparsity
   * pattern used for this
   * matrix.
   */
  std_cxx11::shared_ptr<PointerMatrixBase<VECTOR> > matrix;

  /**
   * Sorted list of pairs denoting
   * the index of the variable and
   * the value to which it shall be
   * fixed.
   */
  std::vector<IndexValuePair> constraints;

  /**
   * Do the pre-filtering step,
   * i.e. zero out those components
   * that belong to constrained
   * degrees of freedom.
   */
  void pre_filter (VECTOR &v) const;

  /**
   * Do the postfiltering step,
   * i.e. set constrained degrees
   * of freedom to the value of the
   * input vector, as the matrix
   * contains only ones on the
   * diagonal for these degrees of
   * freedom.
   */
  void post_filter (const VECTOR &in,
                    VECTOR       &out) const;

  friend class Accessor;
  /**
   * FilteredMatrixBlock accesses
   * pre_filter() and post_filter().
   */
  friend class FilteredMatrixBlock<VECTOR>;
};

/*@}*/
/*---------------------- Inline functions -----------------------------------*/


//--------------------------------Iterators--------------------------------------//

template<class VECTOR>
inline
FilteredMatrix<VECTOR>::Accessor::Accessor(
  const FilteredMatrix<VECTOR> *matrix,
  const size_type index)
  :
  matrix(matrix),
  index(index)
{
  Assert (index <= matrix->constraints.size(),
          ExcIndexRange(index, 0, matrix->constraints.size()));
}



template<class VECTOR>
inline
types::global_dof_index
FilteredMatrix<VECTOR>::Accessor::row() const
{
  return matrix->constraints[index].first;
}



template<class VECTOR>
inline
types::global_dof_index
FilteredMatrix<VECTOR>::Accessor::column() const
{
  return matrix->constraints[index].first;
}



template<class VECTOR>
inline
double
FilteredMatrix<VECTOR>::Accessor::value() const
{
  return matrix->constraints[index].second;
}



template<class VECTOR>
inline
void
FilteredMatrix<VECTOR>::Accessor::advance()
{
  Assert (index < matrix->constraints.size(), ExcIteratorPastEnd());
  ++index;
}




template<class VECTOR>
inline
FilteredMatrix<VECTOR>::const_iterator::const_iterator(
  const FilteredMatrix<VECTOR> *matrix,
  const size_type index)
  :
  accessor(matrix, index)
{}



template<class VECTOR>
inline
typename FilteredMatrix<VECTOR>::const_iterator &
FilteredMatrix<VECTOR>::const_iterator::operator++ ()
{
  accessor.advance();
  return *this;
}


template <typename number>
inline
const typename FilteredMatrix<number>::Accessor &
FilteredMatrix<number>::const_iterator::operator* () const
{
  return accessor;
}


template <typename number>
inline
const typename FilteredMatrix<number>::Accessor *
FilteredMatrix<number>::const_iterator::operator-> () const
{
  return &accessor;
}


template <typename number>
inline
bool
FilteredMatrix<number>::const_iterator::
operator == (const const_iterator &other) const
{
  return (accessor.index == other.accessor.index
          && accessor.matrix == other.accessor.matrix);
}


template <typename number>
inline
bool
FilteredMatrix<number>::const_iterator::
operator != (const const_iterator &other) const
{
  return ! (*this == other);
}



//------------------------------- FilteredMatrix ---------------------------------------//

template <typename number>
inline
typename FilteredMatrix<number>::const_iterator
FilteredMatrix<number>::begin () const
{
  return const_iterator(this, 0);
}


template <typename number>
inline
typename FilteredMatrix<number>::const_iterator
FilteredMatrix<number>::end () const
{
  return const_iterator(this, constraints.size());
}


template <class VECTOR>
inline
bool
FilteredMatrix<VECTOR>::PairComparison::
operator () (const IndexValuePair &i1,
             const IndexValuePair &i2) const
{
  return (i1.first < i2.first);
}



template <class VECTOR>
template <class MATRIX>
inline
void
FilteredMatrix<VECTOR>::initialize (const MATRIX &m, bool ecs)
{
  matrix.reset (new_pointer_matrix_base(m, VECTOR()));

  expect_constrained_source = ecs;
}



template <class VECTOR>
inline
FilteredMatrix<VECTOR>::FilteredMatrix ()
{}



template <class VECTOR>
inline
FilteredMatrix<VECTOR>::FilteredMatrix (const FilteredMatrix &fm)
  :
  Subscriptor(),
  expect_constrained_source(fm.expect_constrained_source),
  matrix(fm.matrix),
  constraints (fm.constraints)
{}



template <class VECTOR>
template <class MATRIX>
inline
FilteredMatrix<VECTOR>::
FilteredMatrix (const MATRIX &m, bool ecs)
{
  initialize (m, ecs);
}



template <class VECTOR>
inline
FilteredMatrix<VECTOR> &
FilteredMatrix<VECTOR>::operator = (const FilteredMatrix &fm)
{
  matrix = fm.matrix;
  expect_constrained_source = fm.expect_constrained_source;
  constraints = fm.constraints;
  return *this;
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::add_constraint (const size_type index, const double value)
{
  // add new constraint to end
  constraints.push_back(IndexValuePair(index, value));
}



template <class VECTOR>
template <class ConstraintList>
inline
void
FilteredMatrix<VECTOR>::add_constraints (const ConstraintList &new_constraints)
{
  // add new constraints to end
  const size_type old_size = constraints.size();
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
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::clear_constraints ()
{
  // swap vectors to release memory
  std::vector<IndexValuePair> empty;
  constraints.swap (empty);
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::clear ()
{
  clear_constraints();
  matrix.reset();
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::apply_constraints (
  VECTOR     &v,
  const bool  /* matrix_is_symmetric */) const
{
  GrowingVectorMemory<VECTOR> mem;
  typename VectorMemory<VECTOR>::Pointer tmp_vector(mem);
  tmp_vector->reinit(v);
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    {
      Assert(numbers::is_finite(i->second), ExcNumberNotFinite());
      (*tmp_vector)(i->first) = -i->second;
    }

  // This vmult is without bc, to get
  // the rhs correction in a correct
  // way.
  matrix->vmult_add(v, *tmp_vector);
  // finally set constrained
  // entries themselves
  for (i=constraints.begin(); i!=e; ++i)
    {
      Assert(numbers::is_finite(i->second), ExcNumberNotFinite());
      v(i->first) = i->second;
    }
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::pre_filter (VECTOR &v) const
{
  // iterate over all constraints and
  // zero out value
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    v(i->first) = 0;
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::post_filter (const VECTOR &in,
                                     VECTOR       &out) const
{
  // iterate over all constraints and
  // set value correctly
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    {
      Assert(numbers::is_finite(in(i->first)), ExcNumberNotFinite());
      out(i->first) = in(i->first);
    }
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::vmult (VECTOR &dst, const VECTOR &src) const
{
  if (!expect_constrained_source)
    {
      GrowingVectorMemory<VECTOR> mem;
      VECTOR *tmp_vector = mem.alloc();
      // first copy over src vector and
      // pre-filter
      tmp_vector->reinit(src, true);
      *tmp_vector = src;
      pre_filter (*tmp_vector);
      // then let matrix do its work
      matrix->vmult (dst, *tmp_vector);
      mem.free(tmp_vector);
    }
  else
    {
      matrix->vmult (dst, src);
    }

  // finally do post-filtering
  post_filter (src, dst);
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  if (!expect_constrained_source)
    {
      GrowingVectorMemory<VECTOR> mem;
      VECTOR *tmp_vector = mem.alloc();
      // first copy over src vector and
      // pre-filter
      tmp_vector->reinit(src, true);
      *tmp_vector = src;
      pre_filter (*tmp_vector);
      // then let matrix do its work
      matrix->Tvmult (dst, *tmp_vector);
      mem.free(tmp_vector);
    }
  else
    {
      matrix->Tvmult (dst, src);
    }

  // finally do post-filtering
  post_filter (src, dst);
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  if (!expect_constrained_source)
    {
      GrowingVectorMemory<VECTOR> mem;
      VECTOR *tmp_vector = mem.alloc();
      // first copy over src vector and
      // pre-filter
      tmp_vector->reinit(src, true);
      *tmp_vector = src;
      pre_filter (*tmp_vector);
      // then let matrix do its work
      matrix->vmult_add (dst, *tmp_vector);
      mem.free(tmp_vector);
    }
  else
    {
      matrix->vmult_add (dst, src);
    }

  // finally do post-filtering
  post_filter (src, dst);
}



template <class VECTOR>
inline
void
FilteredMatrix<VECTOR>::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  if (!expect_constrained_source)
    {
      GrowingVectorMemory<VECTOR> mem;
      VECTOR *tmp_vector = mem.alloc();
      // first copy over src vector and
      // pre-filter
      tmp_vector->reinit(src, true);
      *tmp_vector = src;
      pre_filter (*tmp_vector);
      // then let matrix do its work
      matrix->Tvmult_add (dst, *tmp_vector);
      mem.free(tmp_vector);
    }
  else
    {
      matrix->Tvmult_add (dst, src);
    }

  // finally do post-filtering
  post_filter (src, dst);
}



template <class VECTOR>
inline
std::size_t
FilteredMatrix<VECTOR>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (matrix) +
          MemoryConsumption::memory_consumption (constraints));
}



DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   filtered_matrix.h     ---------------------------*/
