// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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


#ifndef dealii__constraint_matrix_h
#define dealii__constraint_matrix_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/lac/vector.h>

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <complex>


DEAL_II_NAMESPACE_OPEN

template<int dim, class T> class Table;
template <typename> class FullMatrix;
class SparsityPattern;
class DynamicSparsityPattern;
class BlockSparsityPattern;
class BlockDynamicSparsityPattern;
template <typename number> class SparseMatrix;
template <typename number> class BlockSparseMatrix;

namespace internals
{
  class GlobalRowsFromLocal;
}


//TODO[WB]: We should have a function of the kind
//   ConstraintMatrix::add_constraint (const size_type constrained_dof,
//                const std::vector<std::pair<size_type, double> > &entries,
//                const double inhomogeneity = 0);
// rather than building up constraints piecemeal through add_line/add_entry
// etc. This would also eliminate the possibility of accidentally changing
// existing constraints into something pointless, see the discussion on the
// mailing list on "Tiny bug in interpolate_boundary_values" in Sept. 2010.

/**
 * This class implements dealing with linear (possibly inhomogeneous)
 * constraints on degrees of freedom. The concept and origin of such
 * constraints is extensively described in the
 * @ref constraints
 * module. The class is meant to deal with a limited number of constraints
 * relative to the total number of degrees of freedom, for example a few per
 * cent up to maybe 30 per cent; and with a linear combination of <i>M</i>
 * other degrees of freedom where <i>M</i> is also relatively small (no larger
 * than at most around the average number of entries per row of a linear
 * system). It is <em>not</em> meant to describe full rank linear systems.
 *
 * The algorithms used in the implementation of this class are described in
 * some detail in the
 * @ref hp_paper "hp paper".
 * There is also a significant amount of documentation on how to use this
 * class in the
 * @ref constraints
 * module.
 *
 *
 * <h3>Description of constraints</h3>
 *
 * Each "line" in objects of this class corresponds to one constrained degree
 * of freedom, with the number of the line being <i>i</i>, entered by using
 * add_line() or add_lines(). The entries in this line are pairs of the form
 * (<i>j</i>,<i>a<sub>ij</sub></i>), which are added by add_entry() or
 * add_entries(). The organization is essentially a SparsityPattern, but with
 * only a few lines containing nonzero elements, and  therefore no data wasted
 * on the others. For each line, which has been added by the mechanism above,
 * an elimination of the constrained degree of freedom of the form
 * @f[
 *  x_i = \sum_j a_{ij} x_j + b_i
 * @f]
 * is performed, where <i>b<sub>i</sub></i> is optional and set by
 * set_inhomogeneity(). Thus, if a constraint is formulated for instance as a
 * zero mean value of several degrees of freedom, one of the degrees has to be
 * chosen to be eliminated.
 *
 * Note that the constraints are linear in the <i>x<sub>i</sub></i>, and that
 * there might be a constant (non-homogeneous) term in the constraint. This is
 * exactly the form we need for hanging node constraints, where we need to
 * constrain one degree of freedom in terms of others. There are other
 * conditions of this form possible, for example for implementing mean value
 * conditions as is done in the step-11 tutorial program. The name of the
 * class stems from the fact that these constraints can be represented in
 * matrix form as <b>X</b> <i>x</i> = <i>b</i>, and this object then describes
 * the matrix <b>X</b> (and the vector <i>b</i>; originally, the
 * ConstraintMatrix class was only meant to handle homogenous constraints
 * where <i>b</i>=0, thus the name). The most frequent way to create/fill
 * objects of this type is using the DoFTools::make_hanging_node_constraints()
 * function. The use of these objects is first explained in step-6.
 *
 * Objects of the present type are organized in lines (rows), but only those
 * lines are stored where constraints are present. New constraints are added
 * by adding new lines using the add_line() function, and then populating it
 * using the add_entry() function to a given line, or add_entries() to add
 * more than one entry at a time. The right hand side element, if nonzero, can
 * be set using the set_inhomogeneity() function. After all constraints have
 * been added, you need to call close(), which compresses the storage format
 * and sorts the entries.
 *
 * @note Many of the algorithms this class implements are discussed in the
 * @ref hp_paper.
 * The algorithms are also related to those shown in <i>M. S. Shephard: Linear
 * multipoint constraints applied via transformation as part of a direct
 * stiffness assembly process. Int. J. Numer. Meth. Engrg., vol. 20 (1984),
 * pp. 2107-2112.</i>, with the difference that the algorithms shown there
 * completely eliminated constrained degrees of freedom, whereas we usually
 * keep them as part of the linear system.
 *
 * @ingroup dofs
 * @ingroup constraints
 * @author Wolfgang Bangerth, Martin Kronbichler, 1998, 2004, 2008, 2009
 */
class ConstraintMatrix : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * An enum that describes what should happen if the two ConstraintMatrix
   * objects involved in a call to the merge() function happen to have
   * constraints on the same degrees of freedom.
   */
  enum MergeConflictBehavior
  {
    /**
     * Throw an exception if the two objects concerned have conflicting
     * constraints on the same degree of freedom.
     */
    no_conflicts_allowed,

    /**
     * In an operation <code>cm1.merge(cm2)</code>, if <code>cm1</code> and
     * <code>cm2</code> have constraints on the same degree of freedom, take
     * the one from <code>cm1</code>.
     */
    left_object_wins,

    /**
     * In an operation <code>cm1.merge(cm2)</code>, if <code>cm1</code> and
     * <code>cm2</code> have constraints on the same degree of freedom, take
     * the one from <code>cm2</code>.
     */
    right_object_wins
  };

  /**
   * Constructor. The supplied IndexSet defines which indices might be
   * constrained inside this ConstraintMatrix. In a calculation with a
   * DoFHandler object based on parallel::distributed::Triangulation or
   * parallel::shared::Triangulation, one should use the set of locally
   * relevant dofs (see
   * @ref GlossLocallyRelevantDof).
   *
   * The given IndexSet allows the ConstraintMatrix to save memory by just not
   * caring about degrees of freedom that are not of importance to the current
   * processor. Alternatively, if no such IndexSet is provided, internal data
   * structures for <i>all</i> possible indices will be created, leading to
   * memory consumption on every processor that is proportional to the
   * <i>overall</i> size of the problem, not just proportional to the size of
   * the portion of the overall problem that is handled by the current
   * processor.
   */
  explicit ConstraintMatrix (const IndexSet &local_constraints = IndexSet());

  /**
   * Copy constructor
   */
  explicit ConstraintMatrix (const ConstraintMatrix &constraint_matrix);

  /**
   * clear() the ConstraintMatrix object and supply an IndexSet with lines
   * that may be constrained. This function is only relevant in the
   * distributed case to supply a different IndexSet. Otherwise this routine
   * is equivalent to calling clear(). See the constructor for details.
   */
  void reinit (const IndexSet &local_constraints = IndexSet());

  /**
   * Determines if we can store a constraint for the given @p line_index. This
   * routine only matters in the distributed case and checks if the IndexSet
   * allows storage of this line. Always returns true if not in the
   * distributed case.
   */
  bool can_store_line (const size_type line_index) const;

  /**
   * Returns the index set describing locally relevant lines if any are
   * present. Note that if no local lines were given, this represents an empty
   * IndexSet, whereas otherwise it contains the global problem size and the
   * local range.
   */
  const IndexSet &get_local_lines() const;

  /**
   * This function copies the content of @p constraints_in with DoFs that are
   * element of the IndexSet @p filter. Elements that are not present in the
   * IndexSet are ignored. All DoFs will be transformed to local index space
   * of the filter, both the constrained DoFs and the other DoFs these entries
   * are constrained to. The local index space of the filter is a contiguous
   * numbering of all (global) DoFs that are elements in the filter.
   *
   * If, for example, the filter represents the range <tt>[10,20)</tt>, and
   * the constraint matrix @p constraints_in includes the global indices
   * <tt>{7,13,14}</tt>, the indices <tt>{3,4}</tt> are added to the calling
   * constraint matrix (since 13 and 14 are elements in the filter and element
   * 13 is the fourth element in the index, and 14 is the fifth).
   *
   * This function provides an easy way to create a ConstraintMatrix for
   * certain vector components in a vector-valued problem from a full
   * ConstraintMatrix, i.e. extracting a diagonal subblock from a larger
   * ConstraintMatrix. The block is specified by the IndexSet argument.
   */
  void add_selected_constraints (const ConstraintMatrix &constraints_in,
                                 const IndexSet         &filter);

  /**
   * @name Adding constraints
   * @{
   */

  /**
   * Add a new line to the matrix. If the line already exists, then the
   * function simply returns without doing anything.
   */
  void add_line (const size_type line);

  /**
   * Call the first add_line() function for every index <code>i</code> for
   * which <code>lines[i]</code> is true.
   *
   * This function essentially exists to allow adding several constraints of
   * the form <i>x<sub>i</sub></i>=0 all at once, where the set of indices
   * <i>i</i> for which these constraints should be added are given by the
   * argument of this function. On the other hand, just as if the single-
   * argument add_line() function were called repeatedly, the constraints can
   * later be modified to include linear dependencies using the add_entry()
   * function as well as inhomogeneities using set_inhomogeneity().
   */
  void add_lines (const std::vector<bool> &lines);

  /**
   * Call the first add_line() function for every index <code>i</code> that
   * appears in the argument.
   *
   * This function essentially exists to allow adding several constraints of
   * the form <i>x<sub>i</sub></i>=0 all at once, where the set of indices
   * <i>i</i> for which these constraints should be added are given by the
   * argument of this function. On the other hand, just as if the single-
   * argument add_line() function were called repeatedly, the constraints can
   * later be modified to include linear dependencies using the add_entry()
   * function as well as inhomogeneities using set_inhomogeneity().
   */
  void add_lines (const std::set<size_type> &lines);

  /**
   * Call the first add_line() function for every index <code>i</code> that
   * appears in the argument.
   *
   * This function essentially exists to allow adding several constraints of
   * the form <i>x<sub>i</sub></i>=0 all at once, where the set of indices
   * <i>i</i> for which these constraints should be added are given by the
   * argument of this function. On the other hand, just as if the single-
   * argument add_line() function were called repeatedly, the constraints can
   * later be modified to include linear dependencies using the add_entry()
   * function as well as inhomogeneities using set_inhomogeneity().
   */
  void add_lines (const IndexSet &lines);

  /**
   * Add an entry to a given line. The list of lines is searched from the back
   * to the front, so clever programming would add a new line (which is pushed
   * to the back) and immediately afterwards fill the entries of that line.
   * This way, no expensive searching is needed.
   *
   * If an entry with the same indices as the one this function call denotes
   * already exists, then this function simply returns provided that the value
   * of the entry is the same. Thus, it does no harm to enter a constraint
   * twice.
   */
  void add_entry (const size_type line,
                  const size_type column,
                  const double value);

  /**
   * Add a whole series of entries, denoted by pairs of column indices and
   * values, to a line of constraints. This function is equivalent to calling
   * the preceding function several times, but is faster.
   */
  void add_entries (const size_type                                  line,
                    const std::vector<std::pair<size_type,double> > &col_val_pairs);

  /**
   * Set an inhomogeneity to the constraint line <i>i</i>, according to the
   * discussion in the general class description.
   *
   * @note the line needs to be added with one of the add_line() calls first.
   */
  void set_inhomogeneity (const size_type line,
                          const double    value);

  /**
   * Close the filling of entries. Since the lines of a matrix of this type
   * are usually filled in an arbitrary order and since we do not want to use
   * associative constrainers to store the lines, we need to sort the lines
   * and within the lines the columns before usage of the matrix. This is done
   * through this function.
   *
   * Also, zero entries are discarded, since they are not needed.
   *
   * After closing, no more entries are accepted. If the object was already
   * closed, then this function returns immediately.
   *
   * This function also resolves chains of constraints. For example, degree of
   * freedom 13 may be constrained to $u_{13} = \frac{u_3}{2} + \frac{u_7}{2}$
   * while degree of freedom 7 is itself constrained as $u_{7} = \frac{u_2}{2}
   * + \frac{u_4}{2}$. Then, the resolution will be that $u_{13} =
   * \frac{u_3}{2} + \frac{u_2}{4} + \frac{u_4}{4}$. Note, however, that
   * cycles in this graph of constraints are not allowed, i.e. for example
   * $u_4$ may not be constrained, directly or indirectly, to $u_{13}$ again.
   */
  void close ();

  /**
   * Merge the constraints represented by the object given as argument into
   * the constraints represented by this object. Both objects may or may not
   * be closed (by having their function close() called before). If this
   * object was closed before, then it will be closed afterwards as well.
   * Note, however, that if the other argument is closed, then merging may be
   * significantly faster.
   *
   * Using the default value of the second arguments, the constraints in each
   * of the two objects (the old one represented by this object and the
   * argument) may not refer to the same degree of freedom, i.e. a degree of
   * freedom that is constrained in one object may not be constrained in the
   * second. If this is nevertheless the case, an exception is thrown.
   * However, this behavior can be changed by providing a different value for
   * the second argument.
   */
  void merge (const ConstraintMatrix &other_constraints,
              const MergeConflictBehavior merge_conflict_behavior = no_conflicts_allowed);

  /**
   * Shift all entries of this matrix down @p offset rows and over @p offset
   * columns.
   *
   * This function is useful if you are building block matrices, where all
   * blocks are built by the same DoFHandler object, i.e. the matrix size is
   * larger than the number of degrees of freedom. Since several matrix rows
   * and columns correspond to the same degrees of freedom, you'd generate
   * several constraint objects, then shift them, and finally merge() them
   * together again.
   */
  void shift (const size_type offset);

  /**
   * Clear all entries of this matrix. Reset the flag determining whether new
   * entries are accepted or not.
   *
   * This function may be called also on objects which are empty or already
   * cleared.
   */
  void clear ();

  /**
   * @}
   */


  /**
   * @name Querying constraints
   * @{
   */

  /**
   * Return number of constraints stored in this matrix.
   */
  size_type n_constraints () const;

  /**
   * Return whether the degree of freedom with number @p index is a
   * constrained one.
   *
   * Note that if close() was called before, then this function is
   * significantly faster, since then the constrained degrees of freedom are
   * sorted and we can do a binary search, while before close() was called, we
   * have to perform a linear search through all entries.
   */
  bool is_constrained (const size_type index) const;

  /**
   * Return whether the dof is constrained, and whether it is constrained to
   * only one other degree of freedom with weight one. The function therefore
   * returns whether the degree of freedom would simply be eliminated in favor
   * of exactly one other degree of freedom.
   *
   * The function returns @p false if either the degree of freedom is not
   * constrained at all, or if it is constrained to more than one other degree
   * of freedom, or if it is constrained to only one degree of freedom but
   * with a weight different from one.
   */
  bool is_identity_constrained (const size_type index) const;

  /**
   * Return whether the two given degrees of freedom are linked by an equality
   * constraint that either constrains index1 to be so that
   * <code>index1=index2</code> or constrains index2 so that
   * <code>index2=index1</code>.
   */
  bool are_identity_constrained (const size_type index1,
                                 const size_type index2) const;

  /**
   * Return the maximum number of other dofs that one dof is constrained to.
   * For example, in 2d a hanging node is constrained only to its two
   * neighbors, so the returned value would be 2. However, for higher order
   * elements and/or higher dimensions, or other types of constraints, this
   * number is no more obvious.
   *
   * The name indicates that within the system matrix, references to a
   * constrained node are indirected to the nodes it is constrained to.
   */
  size_type max_constraint_indirections () const;

  /**
   * Returns <tt>true</tt> in case the dof is constrained and there is a non-
   * trivial inhomogeneous values set to the dof.
   */
  bool is_inhomogeneously_constrained (const size_type index) const;

  /**
   * Returns <tt>false</tt> if all constraints in the ConstraintMatrix are
   * homogeneous ones, and <tt>true</tt> if there is at least one
   * inhomogeneity.
   */
  bool has_inhomogeneities () const;

  /**
   * Returns a pointer to the the vector of entries if a line is constrained,
   * and a zero pointer in case the dof is not constrained.
   */
  const std::vector<std::pair<size_type,double> > *
  get_constraint_entries (const size_type line) const;

  /**
   * Returns the value of the inhomogeneity stored in the constrained dof @p
   * line. Unconstrained dofs also return a zero value.
   */
  double get_inhomogeneity (const size_type line) const;

  /**
   * Print the constraint lines. Mainly for debugging purposes.
   *
   * This function writes out all entries in the constraint matrix lines with
   * their value in the form <tt>row col : value</tt>. Unconstrained lines
   * containing only one identity entry are not stored in this object and are
   * not printed.
   */
  void print (std::ostream &) const;

  /**
   * Write the graph of constraints in 'dot' format. 'dot' is a program that
   * can take a list of nodes and produce a graphical representation of the
   * graph of constrained degrees of freedom and the degrees of freedom they
   * are constrained to.
   *
   * The output of this function can be used as input to the 'dot' program
   * that can convert the graph into a graphical representation in postscript,
   * png, xfig, and a number of other formats.
   *
   * This function exists mostly for debugging purposes.
   */
  void write_dot (std::ostream &) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  /**
   * Add the constraint indices associated to the indices in the given vector.
   * After a call to this function, the indices vector contains the initial
   * elements and all the associated constrained indices. This function sorts
   * the elements and suppresses duplicates.
   */
  void resolve_indices(std::vector<types::global_dof_index> &indices) const;

  /**
   * @}
   */

  /**
   * @name Eliminating constraints from linear systems after their creation
   * @{
   */


  /**
   * Condense a sparsity pattern. The name of the function mimics the name of
   * the function we use to condense linear systems, but it is a bit of a
   * misnomer for the current context. This is because in the context of
   * linear systems, we eliminate certain rows and columns of the linear
   * system, i.e., we "reduce" or "condense" the linear system. On the other
   * hand, in the current context, the functions does not remove nonzero
   * entries from the sparsity pattern. Rather, it adds those nonzero entry
   * locations to the sparsity pattern that will later be needed for the
   * process of condensation of constrained degrees of freedom from a linear
   * system.
   *
   * Since this function adds new nonzero entries to the sparsity pattern, the
   * given sparsity pattern must not be compressed. The constraint matrix
   * (i.e., the current object) must be closed. The sparsity pattern is
   * compressed at the end of the function.
   */
  void condense (SparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square block sparsity patterns.
   */
  void condense (BlockSparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square compressed sparsity
   * patterns.
   */
  void condense (DynamicSparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square compressed sparsity
   * patterns.
   */
  void condense (BlockDynamicSparsityPattern &sparsity) const;

  /**
   * Condense a given matrix, i.e., eliminate the rows and columns of the
   * matrix that correspond to constrained degrees of freedom.
   *
   * See the general documentation of this class for more detailed
   * information.
   */
  template<typename number>
  void condense (SparseMatrix<number> &matrix) const;

  /**
   * Same function as above, but condenses square block sparse matrices.
   */
  template <typename number>
  void condense (BlockSparseMatrix<number> &matrix) const;

  /**
   * Condense the given vector in-place. The @p VectorType may be a
   * Vector<float>, Vector<double>, BlockVector<tt><...></tt>, a PETSc or
   * Trilinos vector wrapper class, or any other type having the same
   * interface. Note that this function does not take any inhomogeneity into
   * account and throws an exception in case there are any inhomogeneities.
   * Use the function using both a matrix and vector for that case.
   *
   * @note This function does not work for MPI vectors. Use condense() with
   * two vector arguments instead.
   */
  template <class VectorType>
  void condense (VectorType &vec) const;

  /**
   * The function copies and condenses values from @p vec_ghosted into @p
   * output. In a serial code it is equivalent to calling condense (vec). If
   * called in parallel, @p vec_ghosted is supposed to contain ghost elements
   * while @p output should not.
   */
  template <class VectorType>
  void condense (const VectorType &vec_ghosted,
                 VectorType       &output) const;

  /**
   * Condense a given matrix and a given vector by eliminating rows and
   * columns of the linear system that correspond to constrained degrees of
   * freedom. The sparsity pattern associated with the matrix needs to be
   * condensed and compressed.  This function is the appropriate choice for
   * applying inhomogeneous constraints.
   *
   * The constraint matrix object must be closed to call this function.
   *
   * See the general documentation of this class for more detailed
   * information.
   */
  template<typename number, class VectorType>
  void condense (SparseMatrix<number> &matrix,
                 VectorType           &vector) const;

  /**
   * Same function as above, but condenses square block sparse matrices and
   * vectors.
   */
  template <typename number, class BlockVectorType>
  void condense (BlockSparseMatrix<number> &matrix,
                 BlockVectorType           &vector) const;

  /**
   * Sets the values of all constrained DoFs in a vector to zero.  The @p
   * VectorType may be a Vector<float>, Vector<double>,
   * BlockVector<tt><...></tt>, a PETSc or Trilinos vector wrapper class, or
   * any other type having the same interface.
   */
  template <class VectorType>
  void set_zero (VectorType &vec) const;

  /**
   * @}
   */

  /**
   * @name Eliminating constraints from linear systems during their creation
   * @{
   */

  /**
   * This function takes a vector of local contributions (@p local_vector)
   * corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global vector. In most
   * cases, these local contributions will be the result of an integration
   * over a cell or face of a cell. However, as long as @p local_vector and @p
   * local_dof_indices have the same number of elements, this function is
   * happy with whatever it is given.
   *
   * In contrast to the similar function in the DoFAccessor class, this
   * function also takes care of constraints, i.e. if one of the elements of
   * @p local_dof_indices belongs to a constrained node, then rather than
   * writing the corresponding element of @p local_vector into @p
   * global_vector, the element is distributed to the entries in the global
   * vector to which this particular degree of freedom is constrained.
   *
   * Thus, by using this function to distribute local contributions to the
   * global object, one saves the call to the condense function after the
   * vectors and matrices are fully assembled. On the other hand, by
   * consequence, the function does not only write into the entries enumerated
   * by the @p local_dof_indices array, but also (possibly) others as
   * necessary.
   *
   * Note that this function will apply all constraints as if they were
   * homogeneous. For correctly setting inhomogeneous constraints, use the
   * similar function with a matrix argument or the function with both matrix
   * and vector arguments.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global vector allows for
   * simultaneous access and the access is not to rows with the same global
   * index at the same time. This needs to be made sure from the caller's
   * site. There is no locking mechanism inside this method to prevent data
   * races.
   *
   * @param[in] local_vector Vector of local contributions.
   * @param[in] local_dof_indices Local degrees of freedom indices
   * corresponding to the vector of local contributions.
   * @param[out]  global_vector The global vector to which all local
   * contributions will be added.
   */
  template <class InVector, class OutVector>
  void
  distribute_local_to_global (const InVector               &local_vector,
                              const std::vector<size_type> &local_dof_indices,
                              OutVector                    &global_vector) const;

  /**
   * This function takes a vector of local contributions (@p local_vector)
   * corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global vector. In most
   * cases, these local contributions will be the result of an integration
   * over a cell or face of a cell. However, as long as @p local_vector and @p
   * local_dof_indices have the same number of elements, this function is
   * happy with whatever it is given.
   *
   * In contrast to the similar function in the DoFAccessor class, this
   * function also takes care of constraints, i.e. if one of the elements of
   * @p local_dof_indices belongs to a constrained node, then rather than
   * writing the corresponding element of @p local_vector into @p
   * global_vector, the element is distributed to the entries in the global
   * vector to which this particular degree of freedom is constrained.
   *
   * Thus, by using this function to distribute local contributions to the
   * global object, one saves the call to the condense function after the
   * vectors and matrices are fully assembled. On the other hand, by
   * consequence, the function does not only write into the entries enumerated
   * by the @p local_dof_indices array, but also (possibly) others as
   * necessary. This includes writing into diagonal elements of the matrix if
   * the corresponding degree of freedom is constrained.
   *
   * The fourth argument <tt>local_matrix</tt> is intended to be used in case
   * one wants to apply inhomogeneous constraints on the vector only. Such a
   * situation could be where one wants to assemble of a right hand side
   * vector on a problem with inhomogeneous constraints, but the global matrix
   * has been assembled previously. A typical example of this is a time
   * stepping algorithm where the stiffness matrix is assembled once, and the
   * right hand side updated every time step. Note that, however, the entries
   * in the columns of the local matrix have to be exactly the same as those
   * that have been written into the global matrix. Otherwise, this function
   * will not be able to correctly handle inhomogeneities.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global vector allows for
   * simultaneous access and the access is not to rows with the same global
   * index at the same time. This needs to be made sure from the caller's
   * site. There is no locking mechanism inside this method to prevent data
   * races.
   */
  template <typename VectorType, typename LocalType>
  void
  distribute_local_to_global (const Vector<LocalType>      &local_vector,
                              const std::vector<size_type> &local_dof_indices,
                              VectorType                   &global_vector,
                              const FullMatrix<LocalType>  &local_matrix) const;

  /**
   * Same as the previous function, except that it uses two (possibly) different
   * index sets to correctly handle inhomogeneities when the local matrix is
   * computed from a combination of two neighboring elements, for example for an
   * edge integral term in DG. Note that in the case that these two elements have
   * different polynomial degree, the local matrix is rectangular.
   *
   * <tt>local_dof_indices_row</tt> is the set of row indices and
   * <tt>local_dof_indices_col</tt> is the set of column indices of the local matrix.
   * <tt>diagonal=false</tt> says whether the two index sets are equal or not.
   *
   * If both index sets are equal, <tt>diagonal</tt> must be set to true or we
   * simply use the previous function. If both index sets are different (diagonal=false)
   * the <tt>global_vector</tt> is modified to handle inhomogeneities but no
   * entries from <tt>local_vector</tt> are added. Note that the edge integrals for inner
   * edged for DG do not contribute any values to the right hand side.
   */
  template <typename VectorType, typename LocalType>
  void
  distribute_local_to_global (const Vector<LocalType>      &local_vector,
                              const std::vector<size_type> &local_dof_indices_row,
                              const std::vector<size_type> &local_dof_indices_col,
                              VectorType                   &global_vector,
                              const FullMatrix<LocalType>  &local_matrix,
                              bool diagonal = false) const;

  /**
   * Enter a single value into a result vector, obeying constraints.
   */
  template <class VectorType>
  void
  distribute_local_to_global (const size_type index,
                              const double    value,
                              VectorType     &global_vector) const;

  /**
   * This function takes a pointer to a vector of local contributions (@p
   * local_vector) corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global vector. In most
   * cases, these local contributions will be the result of an integration
   * over a cell or face of a cell. However, as long as the entries in @p
   * local_dof_indices indicate reasonable global vector entries, this
   * function is happy with whatever it is given.
   *
   * If one of the elements of @p local_dof_indices belongs to a constrained
   * node, then rather than writing the corresponding element of @p
   * local_vector into @p global_vector, the element is distributed to the
   * entries in the global vector to which this particular degree of freedom
   * is constrained.
   *
   * Thus, by using this function to distribute local contributions to the
   * global object, one saves the call to the condense function after the
   * vectors and matrices are fully assembled. Note that this function
   * completely ignores inhomogeneous constraints.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global vector allows for
   * simultaneous access and the access is not to rows with the same global
   * index at the same time. This needs to be made sure from the caller's
   * site. There is no locking mechanism inside this method to prevent data
   * races.
   */
  template <typename ForwardIteratorVec, typename ForwardIteratorInd,
            class VectorType>
  void
  distribute_local_to_global (ForwardIteratorVec local_vector_begin,
                              ForwardIteratorVec local_vector_end,
                              ForwardIteratorInd local_indices_begin,
                              VectorType        &global_vector) const;

  /**
   * This function takes a matrix of local contributions (@p local_matrix)
   * corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global matrix. In most
   * cases, these local contributions will be the result of an integration
   * over a cell or face of a cell. However, as long as @p local_matrix and @p
   * local_dof_indices have the same number of elements, this function is
   * happy with whatever it is given.
   *
   * In contrast to the similar function in the DoFAccessor class, this
   * function also takes care of constraints, i.e. if one of the elements of
   * @p local_dof_indices belongs to a constrained node, then rather than
   * writing the corresponding element of @p local_matrix into @p
   * global_matrix, the element is distributed to the entries in the global
   * matrix to which this particular degree of freedom is constrained.
   *
   * With this scheme, we never write into rows or columns of constrained
   * degrees of freedom. In order to make sure that the resulting matrix can
   * still be inverted, we need to do something with the diagonal elements
   * corresponding to constrained nodes. Thus, if a degree of freedom in @p
   * local_dof_indices is constrained, we distribute the corresponding entries
   * in the matrix, but also add the absolute value of the diagonal entry of
   * the local matrix to the corresponding entry in the global matrix.
   * Assuming the discretized operator is positive definite, this guarantees
   * that the diagonal entry is always non-zero, positive, and of the same
   * order of magnitude as the other entries of the matrix. On the other hand,
   * when solving a source problem $Au=f$ the exact value of the diagonal
   * element is not important, since the value of the respective degree of
   * freedom will be overwritten by the distribute() call later on anyway.
   *
   * @note The procedure described above adds an unforeseeable number of
   * artificial eigenvalues to the spectrum of the matrix. Therefore, it is
   * recommended to use the equivalent function with two local index vectors
   * in such a case.
   *
   * By using this function to distribute local contributions to the global
   * object, one saves the call to the condense function after the vectors and
   * matrices are fully assembled.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global matrix allows for
   * simultaneous access and the access is not to rows with the same global
   * index at the same time. This needs to be made sure from the caller's
   * site. There is no locking mechanism inside this method to prevent data
   * races.
   */
  template <typename MatrixType>
  void
  distribute_local_to_global (const FullMatrix<typename MatrixType::value_type> &local_matrix,
                              const std::vector<size_type> &local_dof_indices,
                              MatrixType                   &global_matrix) const;

  /**
   * Does almost the same as the function above but can treat general
   * rectangular matrices.  The main difference to achieve this is that the
   * diagonal entries in constrained rows are left untouched instead of being
   * filled with arbitrary values.
   *
   * Since the diagonal entries corresponding to eliminated degrees of freedom
   * are not set, the result may have a zero eigenvalue, if applied to a
   * square matrix. This has to be considered when solving the resulting
   * problems. For solving a source problem $Au=f$, it is possible to set the
   * diagonal entry after building the matrix by a piece of code of the form
   *
   * @code
   *   for (unsigned int i=0;i<matrix.m();++i)
   *     if (constraints.is_constrained(i))
   *       matrix.diag_element(i) = 1.;
   * @endcode
   *
   * The value of one which is used here is arbitrary, but in the context of
   * Krylov space methods uncritical, since it corresponds to an invariant
   * subspace. If the other matrix entries are smaller or larger by a factor
   * close to machine accuracy, it may be advisable to adjust it.
   *
   * For solving eigenvalue problems, this will only add one spurious zero
   * eigenvalue (with a multiplicity that is possibly greater than one).
   * Taking this into account, nothing else has to be changed.
   */
  template <typename MatrixType>
  void
  distribute_local_to_global (const FullMatrix<typename MatrixType::value_type> &local_matrix,
                              const std::vector<size_type> &row_indices,
                              const std::vector<size_type> &col_indices,
                              MatrixType                   &global_matrix) const;

  /**
   * This function simultaneously writes elements into matrix and vector,
   * according to the constraints specified by the calling ConstraintMatrix.
   * This function can correctly handle inhomogeneous constraints as well. For
   * the parameter use_inhomogeneities_for_rhs see the documentation in
   * @ref constraints
   * module.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global matrix and vector allow
   * for simultaneous access and the access is not to rows with the same
   * global index at the same time. This needs to be made sure from the
   * caller's site. There is no locking mechanism inside this method to
   * prevent data races.
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global (const FullMatrix<typename MatrixType::value_type> &local_matrix,
                              const Vector<typename VectorType::value_type>     &local_vector,
                              const std::vector<size_type>  &local_dof_indices,
                              MatrixType                    &global_matrix,
                              VectorType                    &global_vector,
                              bool                          use_inhomogeneities_for_rhs = false) const;

  /**
   * Do a similar operation as the distribute_local_to_global() function that
   * distributes writing entries into a matrix for constrained degrees of
   * freedom, except that here we don't write into a matrix but only allocate
   * sparsity pattern entries.
   *
   * As explained in the
   * @ref hp_paper "hp paper"
   * and in step-27, first allocating a sparsity pattern and later coming back
   * and allocating additional entries for those matrix entries that will be
   * written to due to the elimination of constrained degrees of freedom
   * (using ConstraintMatrix::condense() ), can be a very expensive procedure.
   * It is cheaper to allocate these entries right away without having to do a
   * second pass over the sparsity pattern object. This function does exactly
   * that.
   *
   * Because the function only allocates entries in a sparsity pattern, all it
   * needs to know are the degrees of freedom that couple to each other.
   * Unlike the previous function, no actual values are written, so the second
   * input argument is not necessary here.
   *
   * The third argument to this function, keep_constrained_entries determines
   * whether the function shall allocate entries in the sparsity pattern at
   * all for entries that will later be set to zero upon condensation of the
   * matrix. These entries are necessary if the matrix is built unconstrained,
   * and only later condensed. They are not necessary if the matrix is built
   * using the distribute_local_to_global() function of this class which
   * distributes entries right away when copying a local matrix into a global
   * object. The default of this argument is true, meaning to allocate the few
   * entries that may later be set to zero.
   *
   * By default, the function adds entries for all pairs of indices given in
   * the first argument to the sparsity pattern (unless
   * keep_constrained_entries is false). However, sometimes one would like to
   * only add a subset of all of these pairs. In that case, the last argument
   * can be used which specifies a boolean mask which of the pairs of indices
   * should be considered. If the mask is false for a pair of indices, then no
   * entry will be added to the sparsity pattern for this pair, irrespective
   * of whether one or both of the indices correspond to constrained degrees
   * of freedom.
   *
   * This function is not typically called from user code, but is used in the
   * DoFTools::make_sparsity_pattern() function when passed a constraint
   * matrix object.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global sparsity pattern allows
   * for simultaneous access and the access is not to rows with the same
   * global index at the same time. This needs to be made sure from the
   * caller's site. There is no locking mechanism inside this method to
   * prevent data races.
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                               SparsityPatternType          &sparsity_pattern,
                               const bool                    keep_constrained_entries = true,
                               const Table<2,bool>          &dof_mask                 = default_empty_table) const;

  /**
   * Similar to the other function, but for non-quadratic sparsity patterns.
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global (const std::vector<size_type> &row_indices,
                               const std::vector<size_type> &col_indices,
                               SparsityPatternType          &sparsity_pattern,
                               const bool                    keep_constrained_entries = true,
                               const Table<2,bool>          &dof_mask                 = default_empty_table) const;

  /**
   * This function imports values from a global vector (@p global_vector) by
   * applying the constraints to a vector of local values, expressed in
   * iterator format.  In most cases, the local values will be identified by
   * the local dof values on a cell. However, as long as the entries in @p
   * local_dof_indices indicate reasonable global vector entries, this
   * function is happy with whatever it is given.
   *
   * If one of the elements of @p local_dof_indices belongs to a constrained
   * node, then rather than writing the corresponding element of @p
   * global_vector into @p local_vector, the constraints are resolved as the
   * respective distribute function does, i.e., the local entry is constructed
   * from the global entries to which this particular degree of freedom is
   * constrained.
   *
   * In contrast to the similar function get_dof_values in the DoFAccessor
   * class, this function does not need the constrained values to be correctly
   * set (i.e., distribute to be called).
   */
  template <typename ForwardIteratorVec, typename ForwardIteratorInd,
            class VectorType>
  void
  get_dof_values (const VectorType  &global_vector,
                  ForwardIteratorInd local_indices_begin,
                  ForwardIteratorVec local_vector_begin,
                  ForwardIteratorVec local_vector_end) const;

  /**
   * @}
   */

  /**
   * @name Dealing with constraints after solving a linear system
   * @{
   */

  /**
   * Given a vector, set all constrained degrees of freedom to values so that
   * the constraints are satisfied. For example, if the current object stores
   * the constraint $x_3=\frac 12 x_1 + \frac 12 x_2$, then this function will
   * read the values of $x_1$ and $x_1$ from the given vector and set the
   * element $x_3$ according to this constraints. Similarly, if the current
   * object stores the constraint $x_{42}=208$, then this function will set
   * the 42nd element of the given vector to 208.
   *
   * @note If this function is called with a parallel vector @p vec, then the
   * vector must not contain ghost elements.
   */
  template <class VectorType>
  void distribute (VectorType &vec) const;

  /**
   * @}
   */

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcMatrixIsClosed);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcMatrixNotClosed);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcLineInexistant,
                  size_type,
                  << "The specified line " << arg1
                  << " does not exist.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException4 (ExcEntryAlreadyExists,
                  size_type, size_type, double, double,
                  << "The entry for the indices " << arg1 << " and "
                  << arg2 << " already exists, but the values "
                  << arg3 << " (old) and " << arg4 << " (new) differ "
                  << "by " << (arg4-arg3) << ".");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2 (ExcDoFConstrainedToConstrainedDoF,
                  int, int,
                  << "You tried to constrain DoF " << arg1
                  << " to DoF " << arg2
                  << ", but that one is also constrained. This is not allowed!");
  /**
   * Exception.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcDoFIsConstrainedFromBothObjects,
                  size_type,
                  << "Degree of freedom " << arg1
                  << " is constrained from both object in a merge operation.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcDoFIsConstrainedToConstrainedDoF,
                  size_type,
                  << "In the given argument a degree of freedom is constrained "
                  << "to another DoF with number " << arg1
                  << ", which however is constrained by this object. This is not"
                  << " allowed.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcRowNotStoredHere,
                  size_type,
                  << "The index set given to this constraint matrix indicates "
                  << "constraints for degree of freedom " << arg1
                  << " should not be stored by this object, but a constraint "
                  << "is being added.");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2 (ExcColumnNotStoredHere,
                  size_type,
                  size_type,
                  << "The index set given to this constraint matrix indicates "
                  << "constraints using degree of freedom " << arg2
                  << " should not be stored by this object, but a constraint "
                  << "for degree of freedom " << arg1 <<" uses it.");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2 (ExcIncorrectConstraint,
                  int, int,
                  << "While distributing the constraint for DoF "
                  << arg1 << ", it turns out that one of the processors "
                  << "who own the " << arg2
                  << " degrees of freedom that x_" << arg1
                  << " is constrained against does not know about "
                  << "the constraint on x_" << arg1
                  << ". Did you not initialize the ConstraintMatrix "
                  << "with the appropriate locally_relevant set so "
                  << "that every processor who owns a DoF that constrains "
                  << "another DoF also knows about this constraint?");

private:

  /**
   * This class represents one line of a constraint matrix.
   */
  struct ConstraintLine
  {
    /**
     * A data type in which we store the list of entries that make up the
     * homogenous part of a constraint.
     */
    typedef std::vector<std::pair<size_type,double> > Entries;

    /**
     * Number of this line. Since only very few lines are stored, we can not
     * assume a specific order and have to store the line number explicitly.
     */
    size_type line;

    /**
     * Row numbers and values of the entries in this line.
     *
     * For the reason why we use a vector instead of a map and the
     * consequences thereof, the same applies as what is said for
     * ConstraintMatrix::lines.
     */
    Entries entries;

    /**
     * Value of the inhomogeneity.
     */
    double inhomogeneity;

    /**
     * This operator is a bit weird and unintuitive: it compares the line
     * numbers of two lines. We need this to sort the lines; in fact we could
     * do this using a comparison predicate.  However, this way, it is easier,
     * albeit unintuitive since two lines really have no god-given order
     * relation.
     */
    bool operator < (const ConstraintLine &) const;

    /**
     * This operator is likewise weird: it checks whether the line indices of
     * the two operands are equal, irrespective of the fact that the contents
     * of the line may be different.
     */
    bool operator == (const ConstraintLine &) const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t memory_consumption () const;
  };

  /**
   * Store the lines of the matrix.  Entries are usually appended in an
   * arbitrary order and insertion into a vector is done best at the end, so
   * the order is unspecified after all entries are inserted. Sorting of the
   * entries takes place when calling the <tt>close()</tt> function.
   *
   * We could, instead of using a vector, use an associative array, like a map
   * to store the lines. This, however, would mean a much more fragmented heap
   * since it allocates many small objects, and would additionally make usage
   * of this matrix much slower.
   */
  std::vector<ConstraintLine> lines;

  /**
   * A list of size_type that contains the position of the ConstraintLine of a
   * constrained degree of freedom, or numbers::invalid_size_type if the
   * degree of freedom is not constrained. The numbers::invalid_size_type
   * return value returns thus whether there is a constraint line for a given
   * degree of freedom index. Note that this class has no notion of how many
   * degrees of freedom there really are, so if we check whether there is a
   * constraint line for a given degree of freedom, then this vector may
   * actually be shorter than the index of the DoF we check for.
   *
   * This field exists since when adding a new constraint line we have to
   * figure out whether it already exists. Previously, we would simply walk
   * the unsorted list of constraint lines until we either hit the end or
   * found it. This algorithm is O(N) if N is the number of constraints, which
   * makes it O(N^2) when inserting all constraints. For large problems with
   * many constraints, this could easily take 5-10 per cent of the total run
   * time. With this field, we can save this time since we find any constraint
   * in O(1) time or get to know that it a certain degree of freedom is not
   * constrained.
   *
   * To make things worse, traversing the list of existing constraints
   * requires reads from many different places in memory. Thus, in large 3d
   * applications, the add_line() function showed up very prominently in the
   * overall compute time, mainly because it generated a lot of cache misses.
   * This should also be fixed by using the O(1) algorithm to access the
   * fields of this array.
   *
   * The field is useful in a number of other contexts as well, e.g. when one
   * needs random access to the constraints as in all the functions that apply
   * constraints on the fly while add cell contributions into vectors and
   * matrices.
   */
  std::vector<size_type> lines_cache;

  /**
   * This IndexSet is used to limit the lines to save in the ConstraintMatrix
   * to a subset. This is necessary, because the lines_cache vector would
   * become too big in a distributed calculation.
   */
  IndexSet local_lines;

  /**
   * Store whether the arrays are sorted.  If so, no new entries can be added.
   */
  bool sorted;

  /**
   * Internal function to calculate the index of line @p line in the vector
   * lines_cache using local_lines.
   */
  size_type calculate_line_index (const size_type line) const;

  /**
   * Return @p true if the weight of an entry (the second element of the pair)
   * equals zero. This function is used to delete entries with zero weight.
   */
  static bool check_zero_weight (const std::pair<size_type, double> &p);

  /**
   * Dummy table that serves as default argument for function
   * <tt>add_entries_local_to_global()</tt>.
   */
  static const Table<2,bool> default_empty_table;

  /**
   * This function actually implements the local_to_global function for
   * standard (non-block) matrices.
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global (const FullMatrix<typename MatrixType::value_type>  &local_matrix,
                              const Vector<typename VectorType::value_type>      &local_vector,
                              const std::vector<size_type> &local_dof_indices,
                              MatrixType                   &global_matrix,
                              VectorType                   &global_vector,
                              bool                          use_inhomogeneities_for_rhs,
                              internal::bool2type<false>) const;

  /**
   * This function actually implements the local_to_global function for block
   * matrices.
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global (const FullMatrix<typename MatrixType::value_type>  &local_matrix,
                              const Vector<typename VectorType::value_type>      &local_vector,
                              const std::vector<size_type> &local_dof_indices,
                              MatrixType                   &global_matrix,
                              VectorType                   &global_vector,
                              bool                          use_inhomogeneities_for_rhs,
                              internal::bool2type<true>) const;

  /**
   * This function actually implements the local_to_global function for
   * standard (non-block) sparsity types.
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                               SparsityPatternType          &sparsity_pattern,
                               const bool                    keep_constrained_entries,
                               const Table<2,bool>          &dof_mask,
                               internal::bool2type<false>) const;

  /**
   * This function actually implements the local_to_global function for block
   * sparsity types.
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                               SparsityPatternType          &sparsity_pattern,
                               const bool                    keep_constrained_entries,
                               const Table<2,bool>          &dof_mask,
                               internal::bool2type<true>) const;

  /**
   * Internal helper function for distribute_local_to_global function.
   *
   * Creates a list of affected global rows for distribution, including the
   * local rows where the entries come from. The list is sorted according to
   * the global row indices.
   */
  void
  make_sorted_row_list (const std::vector<size_type>   &local_dof_indices,
                        internals::GlobalRowsFromLocal &global_rows) const;

  /**
   * Internal helper function for add_entries_local_to_global function.
   *
   * Creates a list of affected rows for distribution without any additional
   * information, otherwise similar to the other make_sorted_row_list()
   * function.
   */
  void
  make_sorted_row_list (const std::vector<size_type> &local_dof_indices,
                        std::vector<size_type>       &active_dofs) const;

  /**
   * Internal helper function for distribute_local_to_global function.
   */
  template <typename LocalType>
  LocalType
  resolve_vector_entry (const size_type                       i,
                        const internals::GlobalRowsFromLocal &global_rows,
                        const Vector<LocalType>              &local_vector,
                        const std::vector<size_type>         &local_dof_indices,
                        const FullMatrix<LocalType>          &local_matrix) const;

  /**
   * The assignment operator is not implemented for performance reasons. You
   * can clear() or reinit() and merge() manually if needed.
   */
  ConstraintMatrix &operator= (const ConstraintMatrix &other);
};



/* ---------------- template and inline functions ----------------- */

inline
ConstraintMatrix::ConstraintMatrix (const IndexSet &local_constraints)
  :
  lines (),
  local_lines (local_constraints),
  sorted (false)
{
  // make sure the IndexSet is compressed. Otherwise this can lead to crashes
  // that are hard to find (only happen in release mode).
  // see tests/mpi/constraint_matrix_crash_01
  local_lines.compress();
}



inline
ConstraintMatrix::ConstraintMatrix (const ConstraintMatrix &constraint_matrix)
  :
  Subscriptor (),
  lines (constraint_matrix.lines),
  lines_cache (constraint_matrix.lines_cache),
  local_lines (constraint_matrix.local_lines),
  sorted (constraint_matrix.sorted)
{}


inline
void
ConstraintMatrix::add_line (const size_type line)
{
  Assert (sorted==false, ExcMatrixIsClosed());

  // the following can happen when we compute with distributed meshes and dof
  // handlers and we constrain a degree of freedom whose number we don't have
  // locally. if we don't abort here the program will try to allocate several
  // terabytes of memory to resize the various arrays below :-)
  Assert (line != numbers::invalid_size_type,
          ExcInternalError());
  const size_type line_index = calculate_line_index (line);

  // check whether line already exists; it may, in which case we can just quit
  if (is_constrained(line))
    return;

  // if necessary enlarge vector of existing entries for cache
  if (line_index >= lines_cache.size())
    lines_cache.resize (std::max(2*static_cast<size_type>(lines_cache.size()),
                                 line_index+1),
                        numbers::invalid_size_type);

  // push a new line to the end of the list
  lines.push_back (ConstraintLine());
  lines.back().line = line;
  lines.back().inhomogeneity = 0.;
  lines_cache[line_index] = lines.size()-1;
}



inline
void
ConstraintMatrix::add_entry (const size_type line,
                             const size_type column,
                             const double    value)
{
  Assert (sorted==false, ExcMatrixIsClosed());
  Assert (line != column,
          ExcMessage ("Can't constrain a degree of freedom to itself"));

  // if in debug mode, check whether an entry for this column already exists
  // and if it's the same as the one entered at present
  //
  // in any case: exit the function if an entry for this column already
  // exists, since we don't want to enter it twice
  Assert (lines_cache[calculate_line_index(line)] != numbers::invalid_size_type,
          ExcInternalError());
  Assert (!local_lines.size() || local_lines.is_element(column),
          ExcColumnNotStoredHere(line, column));
  ConstraintLine *line_ptr = &lines[lines_cache[calculate_line_index(line)]];
  Assert (line_ptr->line == line, ExcInternalError());
  for (ConstraintLine::Entries::const_iterator
       p=line_ptr->entries.begin();
       p != line_ptr->entries.end(); ++p)
    if (p->first == column)
      {
        Assert (std::fabs(p->second - value) < 1.e-14,
                ExcEntryAlreadyExists(line, column, p->second, value));
        return;
      }

  line_ptr->entries.push_back (std::make_pair(column,value));
}



inline
void
ConstraintMatrix::set_inhomogeneity (const size_type line,
                                     const double    value)
{
  const size_type line_index = calculate_line_index(line);
  Assert( line_index < lines_cache.size() &&
          lines_cache[line_index] != numbers::invalid_size_type,
          ExcMessage("call add_line() before calling set_inhomogeneity()"));
  Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
  ConstraintLine *line_ptr = &lines[lines_cache[line_index]];
  line_ptr->inhomogeneity = value;
}



inline
types::global_dof_index
ConstraintMatrix::n_constraints () const
{
  return lines.size();
}



inline
bool
ConstraintMatrix::is_constrained (const size_type index) const
{
  const size_type line_index = calculate_line_index(index);
  return ((line_index < lines_cache.size())
          &&
          (lines_cache[line_index] != numbers::invalid_size_type));
}



inline
bool
ConstraintMatrix::is_inhomogeneously_constrained (const size_type index) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(index);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return false;
  else
    {
      Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
      return !(lines[lines_cache[line_index]].inhomogeneity == 0);
    }
}



inline
const std::vector<std::pair<types::global_dof_index,double> > *
ConstraintMatrix::get_constraint_entries (const size_type line) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return 0;
  else
    return &lines[lines_cache[line_index]].entries;
}



inline
double
ConstraintMatrix::get_inhomogeneity (const size_type line) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return 0;
  else
    return lines[lines_cache[line_index]].inhomogeneity;
}



inline types::global_dof_index
ConstraintMatrix::calculate_line_index (const size_type line) const
{
  //IndexSet is unused (serial case)
  if (!local_lines.size())
    return line;

  Assert(local_lines.is_element(line),
         ExcRowNotStoredHere(line));

  return local_lines.index_within_set(line);
}



inline bool
ConstraintMatrix::can_store_line (size_type line_index) const
{
  return !local_lines.size() || local_lines.is_element(line_index);
}



inline
const IndexSet &
ConstraintMatrix::get_local_lines () const
{
  return local_lines;
}



template <class VectorType>
inline
void ConstraintMatrix::distribute_local_to_global (
  const size_type index,
  const double    value,
  VectorType     &global_vector) const
{
  Assert (lines.empty() || sorted == true, ExcMatrixNotClosed());

  if (is_constrained(index) == false)
    global_vector(index) += value;
  else
    {
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(index)]];
      for (size_type j=0; j<position.entries.size(); ++j)
        global_vector(position.entries[j].first)
        += value * position.entries[j].second;
    }
}


template <typename ForwardIteratorVec, typename ForwardIteratorInd,
          class VectorType>
inline
void ConstraintMatrix::distribute_local_to_global (
  ForwardIteratorVec local_vector_begin,
  ForwardIteratorVec local_vector_end,
  ForwardIteratorInd local_indices_begin,
  VectorType        &global_vector) const
{
  Assert (lines.empty() || sorted == true, ExcMatrixNotClosed());
  for ( ; local_vector_begin != local_vector_end;
        ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        global_vector(*local_indices_begin) += *local_vector_begin;
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          for (size_type j=0; j<position.entries.size(); ++j)
            global_vector(position.entries[j].first)
            += *local_vector_begin * position.entries[j].second;
        }
    }
}


template <class InVector, class OutVector>
inline
void
ConstraintMatrix::distribute_local_to_global (
  const InVector               &local_vector,
  const std::vector<size_type> &local_dof_indices,
  OutVector                    &global_vector) const
{
  Assert (local_vector.size() == local_dof_indices.size(),
          ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  distribute_local_to_global (local_vector.begin(), local_vector.end(),
                              local_dof_indices.begin(), global_vector);
}



template <typename ForwardIteratorVec, typename ForwardIteratorInd,
          class VectorType>
inline
void ConstraintMatrix::get_dof_values (const VectorType  &global_vector,
                                       ForwardIteratorInd local_indices_begin,
                                       ForwardIteratorVec local_vector_begin,
                                       ForwardIteratorVec local_vector_end) const
{
  Assert (lines.empty() || sorted == true, ExcMatrixNotClosed());
  for ( ; local_vector_begin != local_vector_end;
        ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        *local_vector_begin = global_vector(*local_indices_begin);
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          typename VectorType::value_type value = position.inhomogeneity;
          for (size_type j=0; j<position.entries.size(); ++j)
            value += (global_vector(position.entries[j].first) *
                      position.entries[j].second);
          *local_vector_begin = value;
        }
    }
}


template <typename MatrixType> class BlockMatrixBase;
template <typename SparsityPatternType> class BlockSparsityPatternBase;
template <typename number>     class BlockSparseMatrixEZ;

/**
 * A class that can be used to determine whether a given type is a block
 * matrix type or not. For example,
 * @code
 *   IsBlockMatrix<SparseMatrix<double> >::value
 * @endcode
 * has the value false, whereas
 * @code
 *   IsBlockMatrix<BlockSparseMatrix<double> >::value
 * @endcode
 * is true. This is sometimes useful in template contexts where we may want to
 * do things differently depending on whether a template type denotes a
 * regular or a block matrix type.
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, 2009
 */
template <typename MatrixType>
struct IsBlockMatrix
{
private:
  struct yes_type
  {
    char c[1];
  };
  struct no_type
  {
    char c[2];
  };

  /**
   * Overload returning true if the class is derived from BlockMatrixBase,
   * which is what block matrices do (with the exception of
   * BlockSparseMatrixEZ).
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockMatrixBase<T> *);

  /**
   * Overload returning true if the class is derived from
   * BlockSparsityPatternBase, which is what block sparsity patterns do.
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockSparsityPatternBase<T> *);

  /**
   * Overload for BlockSparseMatrixEZ, which is the only block matrix not
   * derived from BlockMatrixBase at the time of writing this class.
   */
  template <typename T>
  static yes_type check_for_block_matrix (const BlockSparseMatrixEZ<T> *);

  /**
   * Catch all for all other potential matrix types that are not block
   * matrices.
   */
  static no_type check_for_block_matrix (...);

public:
  /**
   * A statically computable value that indicates whether the template
   * argument to this class is a block matrix (in fact whether the type is
   * derived from BlockMatrixBase<T>).
   */
  static const bool value = (sizeof(check_for_block_matrix
                                    ((MatrixType *)0))
                             ==
                             sizeof(yes_type));
};


// instantiation of the static member
template <typename MatrixType>
const bool IsBlockMatrix<MatrixType>::value;


template <typename MatrixType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<typename MatrixType::value_type>     &local_matrix,
                            const std::vector<size_type> &local_dof_indices,
                            MatrixType                   &global_matrix) const
{
  // create a dummy and hand on to the function actually implementing this
  // feature in the cm.templates.h file.
  Vector<typename MatrixType::value_type> dummy(0);
  distribute_local_to_global (local_matrix, dummy, local_dof_indices,
                              global_matrix, dummy, false,
                              dealii::internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}




template <typename MatrixType, typename VectorType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<typename MatrixType::value_type>     &local_matrix,
                            const Vector<typename VectorType::value_type>         &local_vector,
                            const std::vector<size_type> &local_dof_indices,
                            MatrixType                   &global_matrix,
                            VectorType                   &global_vector,
                            bool                          use_inhomogeneities_for_rhs) const
{
  // enter the internal function with the respective block information set,
  // the actual implementation follows in the cm.templates.h file.
  distribute_local_to_global (local_matrix, local_vector, local_dof_indices,
                              global_matrix, global_vector, use_inhomogeneities_for_rhs,
                              dealii::internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}




template <typename SparsityPatternType>
inline
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                             SparsityPatternType          &sparsity_pattern,
                             const bool                    keep_constrained_entries,
                             const Table<2,bool>          &dof_mask) const
{
  // enter the internal function with the respective block information set,
  // the actual implementation follows in the cm.templates.h file.
  add_entries_local_to_global (local_dof_indices, sparsity_pattern,
                               keep_constrained_entries, dof_mask,
                               internal::bool2type<IsBlockMatrix<SparsityPatternType>::value>());
}


DEAL_II_NAMESPACE_CLOSE

#endif
