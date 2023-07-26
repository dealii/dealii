// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_affine_constraints_h
#define dealii_affine_constraints_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/lac/sparsity_pattern_base.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_element_access.h>

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class FullMatrix;
class SparsityPattern;
class DynamicSparsityPattern;
class BlockSparsityPattern;
class BlockDynamicSparsityPattern;
template <typename number>
class SparseMatrix;
template <typename number>
class BlockSparseMatrix;

namespace internal
{
  namespace AffineConstraints
  {
    using size_type = types::global_dof_index;

    /**
     * This struct contains all the information we need to store about each of
     * the global entries (global_row): are they obtained directly by some local
     * entry (local_row) or some constraints (constraint_position). This is not
     * directly used in the user code, but accessed via the GlobalRowsFromLocal.
     *
     * The actions performed here correspond to reshaping the constraint
     * information from global degrees of freedom to local ones (i.e.,
     * cell-related DoFs), and also transforming the constraint information from
     * compressed row storage (each local dof that is constrained has a list of
     * constraint entries associated to it) into compressed column storage based
     * on the cell-related DoFs (we have a list of global degrees of freedom,
     * and to each we have a list of local rows where the entries come from). To
     * increase the speed, we additionally store whether an entry is generated
     * directly from the local degrees of freedom or whether it comes from a
     * constraint.
     */
    struct Distributing
    {
      Distributing(const size_type global_row = numbers::invalid_size_type,
                   const size_type local_row  = numbers::invalid_size_type);

      Distributing(const Distributing &in);

      Distributing &
      operator=(const Distributing &in);

      bool
      operator<(const Distributing &in) const
      {
        return global_row < in.global_row;
      }

      size_type         global_row;
      size_type         local_row;
      mutable size_type constraint_position;
    };



    /**
     * This class represents a cache for constraints that are encountered on a
     * local level. The functionality is similar to
     * std::vector<std::vector<std::pair<uint,double> > >, but tuned so that
     * frequent memory allocation for each entry is avoided. The data is put
     * into a std::vector<std::pair<uint,double> > and the row length is kept
     * fixed at row_length. Both the number of rows and the row length can
     * change is this structure is filled. In that case, the data is
     * rearranged. This is not directly used in the user code, but accessed
     * via the GlobalRowsFromLocal.
     */
    template <typename number>
    struct DataCache
    {
      DataCache();

      void
      reinit();

      size_type
      insert_new_index(const std::pair<size_type, number> &pair);

      void
      append_index(const size_type                     index,
                   const std::pair<size_type, number> &pair);

      size_type
      get_size(const size_type index) const;

      const std::pair<size_type, number> *
      get_entry(const size_type index) const;

      size_type row_length;

      std::vector<std::pair<size_type, number>> data;

      std::vector<size_type> individual_size;
    };



    /**
     * A data structure that collects all the global rows from a local
     * contribution (cell) and their origin (direct/constraint). This
     * is basically a vector consisting of "Distributing" structs
     * using access via the DataCache. The structure provides some
     * specialized sort and insert functions.
     *
     * In case there are no constraints, this is basically a list of pairs
     * `<uint,uint>` with the first index being the global index and the second
     * index the local index. The list is sorted with respect to the global
     * index.
     *
     * In case there are constraints, a global dof might get a contribution also
     * because it gets data from a constrained dof. This means that a global dof
     * might also have indirect contributions from a local dof via a constraint,
     * besides the direct ones.
     *
     * The actions performed here correspond to reshaping the constraint
     * information from global degrees of freedom to local ones (i.e.,
     * cell-related DoFs), and also transforming the constraint information from
     * compressed row storage (each local dof that is constrained has a list of
     * constraint entries associated to it) into compressed column storage based
     * on the cell-related DoFs (we have a list of global degrees of freedom,
     * and to each we have a list of local rows where the entries come from). To
     * increase the speed, we additionally store whether an entry is generated
     * directly from the local degrees of freedom or whether it comes from a
     * constraint.
     */
    template <typename number>
    class GlobalRowsFromLocal
    {
    public:
      /**
       * Constructor.
       */
      GlobalRowsFromLocal();

      void
      reinit(const size_type n_local_rows);

      void
      insert_index(const size_type global_row,
                   const size_type local_row,
                   const number    constraint_value);
      void
      sort();

      void
      print(std::ostream &os);

      /**
       * Return the number of global indices in the struct.
       */
      size_type
      size() const;

      /**
       * Return the number of constraints that are associated to the
       * counter_index-th entry in the list.
       */
      size_type
      size(const size_type counter_index) const;

      /**
       * Return the global row of the counter_index-th entry in the list.
       */
      size_type
      global_row(const size_type counter_index) const;

      /**
       * Return the global row of the counter_index-th entry in the list.
       */
      size_type &
      global_row(const size_type counter_index);

      /**
       * Return the local row in the cell matrix associated with the
       * counter_index-th entry in the list. Return invalid_size_type for
       * constrained rows.
       */
      size_type
      local_row(const size_type counter_index) const;

      /**
       * Return a reference instead of the value as in the function above.
       */
      size_type &
      local_row(const size_type counter_index);

      /**
       * Return the local row in the cell matrix associated with the
       * counter_index-th entry in the list in the index_in_constraint-th
       * position of constraints.
       */
      size_type
      local_row(const size_type counter_index,
                const size_type index_in_constraint) const;

      /**
       * Return the value of the constraint in the counter_index-th entry in
       * the list in the index_in_constraint-th position of constraints.
       */
      number
      constraint_value(const size_type counter_index,
                       const size_type index_in_constraint) const;

      /**
       * Return whether there is one row with indirect contributions (i.e.,
       * there has been at least one constraint with non-trivial
       * ConstraintLine).
       */
      bool
      have_indirect_rows() const;

      /**
       * Append an entry that is constrained. This means that there is one less
       * nontrivial row.
       */
      void
      insert_constraint(const size_type constrained_local_dof);

      /**
       * Return the number of constrained dofs in the structure. Constrained
       * dofs do not contribute directly to the matrix, but are needed in order
       * to set matrix diagonals and resolve inhomogeneities.
       */
      size_type
      n_constraints() const;

      /**
       * Return the number of constrained dofs in the structure that have an
       * inhomogeneity.
       */
      size_type
      n_inhomogeneities() const;

      /**
       * This function tells the structure that the ith constraint is
       * inhomogeneous. inhomogeneous constraints contribute to right hand
       * sides, so to have fast access to them, put them before homogeneous
       * constraints.
       */
      void
      set_ith_constraint_inhomogeneous(const size_type i);

      /**
       * The local row where constraint number i was detected, to find that row
       * easily when the GlobalRowsToLocal has been set up.
       */
      size_type
      constraint_origin(size_type i) const;

      /**
       * A vector that contains all the global ids and the corresponding local
       * ids as well as a pointer to that data where we store how to resolve
       * constraints.
       */
      std::vector<Distributing> total_row_indices;

    private:
      /**
       * A data structure that holds the actual data from the constraints.
       */
      DataCache<number> data_cache;

      /**
       * A number that states how many rows there are, constraints
       * disregarded.
       */
      size_type n_active_rows;

      /**
       * A number that represents the number of rows with
       * inhomogeneous constraints.
       */
      size_type n_inhomogeneous_rows;
    };



    /**
     * Scratch data that is used during calls to distribute_local_to_global and
     * add_entries_local_to_global. In order to avoid frequent memory
     * allocation, we keep the data alive from one call to the next in a static
     * variable. Since we want to allow for different number types in matrices,
     * this is a template.
     *
     * Since each thread gets its private version of scratch data out of a
     * ThreadLocalStorage, no conflicting access can occur. For this to be
     * valid, we need to make sure that no call within
     * distribute_local_to_global is made that by itself can spawn tasks.
     * Otherwise, we might end up in a situation where several threads fight for
     * the data.
     *
     * Access to the scratch data is only through an accessor class which
     * handles the access as well as marks the data as used.
     */
    template <typename number>
    struct ScratchData
    {
      /**
       * Constructor, does nothing.
       */
      ScratchData()
        : in_use(false)
      {}

      /**
       * Copy constructor, does nothing
       */
      ScratchData(const ScratchData &)
        : in_use(false)
      {}

      /**
       * Stores whether the data is currently in use.
       */
      bool in_use;

      /**
       * Temporary array for pairs of indices
       */
      std::vector<std::pair<size_type, size_type>> new_entries;

      /**
       * Temporary array for row indices
       */
      std::vector<size_type> rows;

      /**
       * Temporary array for column indices
       */
      std::vector<size_type> columns;

      /**
       * Temporary array for column values
       */
      std::vector<number> values;

      /**
       * Temporary array for block start indices
       */
      std::vector<size_type> block_starts;

      /**
       * Temporary array for vector indices
       */
      std::vector<size_type> vector_indices;

      /**
       * Temporary array for vector values
       */
      std::vector<number> vector_values;

      /**
       * Data array for reorder row/column indices.
       */
      GlobalRowsFromLocal<number> global_rows;

      /**
       * Data array for reorder row/column indices.
       */
      GlobalRowsFromLocal<number> global_columns;
    };
  } // namespace AffineConstraints
} // namespace internal

namespace internal
{
  namespace AffineConstraintsImplementation
  {
    template <typename VectorType>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 VectorType &                                vec);

    template <class T>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 dealii::Vector<T> &                         vec);

    template <class T>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 dealii::BlockVector<T> &                    vec);
  } // namespace AffineConstraintsImplementation
} // namespace internal
#endif

// TODO[WB]: We should have a function of the kind
//   AffineConstraints::add_constraint (const size_type constrained_dof,
//     const std::vector<std::pair<size_type, number> > &entries,
//     const number inhomogeneity = 0);
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
 * @ref hp_paper "hp-paper".
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
 * the matrix <b>X</b> and the vector <i>b</i>. The most frequent way to
 * create/fill objects of this type is using the
 * DoFTools::make_hanging_node_constraints() function. The use of these
 * objects is first explained in step-6.
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
 */
template <typename number = double>
class AffineConstraints : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * An enum that describes what should happen if the two AffineConstraints
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
   * constrained inside this AffineConstraints container. In a calculation
   * with a DoFHandler object based on parallel::distributed::Triangulation
   * or parallel::shared::Triangulation, one should use the set of locally
   * relevant dofs (see
   * @ref GlossLocallyRelevantDof).
   *
   * The given IndexSet allows the AffineConstraints container to save
   * memory by just not caring about degrees of freedom that are not of
   * importance to the current processor. Alternatively, if no such
   * IndexSet is provided, internal data structures for <i>all</i> possible
   * indices will be created, leading to memory consumption on every
   * processor that is proportional to the <i>overall</i> size of the
   * problem, not just proportional to the size of the portion of the
   * overall problem that is handled by the current processor.
   */
  explicit AffineConstraints(const IndexSet &local_constraints = IndexSet());

  /**
   * Copy constructor
   */
  explicit AffineConstraints(const AffineConstraints &affine_constraints);

  /**
   * Move constructor
   */
  AffineConstraints(AffineConstraints &&affine_constraints) noexcept =
    default; // NOLINT

  /**
   * Copy operator. Like for many other large objects, this operator
   * is deleted to avoid its inadvertent use in places such as
   * accidentally declaring a @p AffineConstraints object as a
   * function argument by value, rather than by reference.
   *
   * However, you can use the copy_from() function to explicitly
   * copy AffineConstraints objects.
   */
  AffineConstraints &
  operator=(const AffineConstraints &) = delete;

  /**
   * Move assignment operator
   */
  AffineConstraints &
  operator=(AffineConstraints &&affine_constraints) noexcept =
    default; // NOLINT

  /**
   * Copy the given object to the current one.
   *
   * This function exists because @p operator=() is explicitly
   * disabled.
   */
  template <typename other_number>
  void
  copy_from(const AffineConstraints<other_number> &other);

  /**
   * clear() the AffineConstraints object and supply an IndexSet with lines
   * that may be constrained. This function is only relevant in the
   * distributed case to supply a different IndexSet. Otherwise this routine
   * is equivalent to calling clear(). See the constructor for details.
   */
  void
  reinit(const IndexSet &local_constraints = IndexSet());

  /**
   * Determines if we can store a constraint for the given @p line_n. This
   * routine only matters in the distributed case and checks if the IndexSet
   * allows storage of this line. Always returns true if not in the
   * distributed case.
   */
  bool
  can_store_line(const size_type line_n) const;

  /**
   * Return the index set describing locally relevant lines if any are
   * present. Note that if no local lines were given, this represents an empty
   * IndexSet, whereas otherwise it contains the global problem size and the
   * local range.
   */
  const IndexSet &
  get_local_lines() const;

  /**
   * This function copies the content of @p constraints_in with DoFs that are
   * element of the IndexSet @p filter. Elements that are not present in the
   * IndexSet are ignored. All DoFs will be transformed to local index space
   * of the filter, both the constrained DoFs and the other DoFs these entries
   * are constrained to. The local index space of the filter is a contiguous
   * numbering of all (global) DoFs that are elements in the filter.
   *
   * If, for example, the filter represents the range <tt>[10,20)</tt>, and
   * the constraints object @p constraints_in includes the global indices
   * <tt>{7,13,14}</tt>, the indices <tt>{3,4}</tt> are added to the calling
   * constraints object (since 13 and 14 are elements in the filter and element
   * 13 is the fourth element in the index, and 14 is the fifth).
   *
   * This function provides an easy way to create a AffineConstraints for
   * certain vector components in a vector-valued problem from a full
   * AffineConstraints, i.e. extracting a diagonal subblock from a larger
   * AffineConstraints. The block is specified by the IndexSet argument.
   */
  void
  add_selected_constraints(const AffineConstraints &constraints_in,
                           const IndexSet &         filter);

  /**
   * @name Adding constraints
   * @{
   */

  /**
   * Add a new line to the matrix. If the line already exists, then the
   * function simply returns without doing anything.
   */
  void
  add_line(const size_type line_n);

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
  void
  add_lines(const std::vector<bool> &lines);

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
  void
  add_lines(const std::set<size_type> &lines);

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
  void
  add_lines(const IndexSet &lines);

  /**
   * Add an entry to a given line. In other words, this function adds
   * a term $a_{ij} x_j$ to the constraints for the $i$th degree of freedom.
   *
   * If an entry with the same indices as the one this function call denotes
   * already exists, then this function simply returns provided that the value
   * of the entry is the same. Thus, it does no harm to enter a constraint
   * twice.
   *
   * @param[in] constrained_dof_index The index $i$ of the degree of freedom
   *   that is being constrained.
   * @param[in] column The index $j$ of the degree of freedom being entered
   *   into the constraint for degree of freedom $i$.
   * @param[in] weight The factor $a_{ij}$ that multiplies $x_j$.
   */
  void
  add_entry(const size_type constrained_dof_index,
            const size_type column,
            const number    weight);

  /**
   * Add a whole series of entries, denoted by pairs of column indices and
   * weight values, to a line of constraints. This function is equivalent to
   * calling the preceding function several times, but is faster.
   */
  void
  add_entries(
    const size_type                                  constrained_dof_index,
    const std::vector<std::pair<size_type, number>> &col_weight_pairs);

  /**
   * Set an inhomogeneity to the constraint for a degree of freedom. In other
   * words, it adds a constant $b_i$ to the constraint for degree of freedom
   * $i$. For this to work, you need to call add_line() first for the given
   * degree of freedom.
   *
   * @param[in] constrained_dof_index The index $i$ of the degree of freedom
   *   that is being constrained.
   * @param[in] value The right hand side value $b_i$ for the constraint on
   *   the degree of freedom $i$.
   */
  void
  set_inhomogeneity(const size_type constrained_dof_index, const number value);

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
   * cycles in this graph of constraints are not allowed, i.e., for example
   * $u_4$ may not itself be constrained, directly or indirectly, to $u_{13}$
   * again.
   */
  void
  close();

  /**
   * Check if the function close() was called or there are no
   * constraints locally, which is normally the case if a dummy
   * AffineConstraints was created for the DG case.
   */
  bool
  is_closed() const;

  /**
   * Check if the function close() was called or there are no
   * constraints globally, which is normally the case if a dummy
   * AffineConstraints was created for the DG case.
   */
  bool
  is_closed(const MPI_Comm comm) const;

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
   *
   * By default, merging two AffineConstraints objects that are initialized
   * with different IndexSet objects is not allowed.
   * This behavior can be altered by setting @p allow_different_local_lines
   * appropriately.
   *
   * Merging a AffineConstraints that is initialized with an IndexSet
   * and one that is not initialized with an IndexSet is not yet implemented.
   */
  template <typename other_number>
  void
  merge(
    const AffineConstraints<other_number> &other_constraints,
    const MergeConflictBehavior merge_conflict_behavior = no_conflicts_allowed,
    const bool                  allow_different_local_lines = false);

  /**
   * Shift all entries of this matrix down @p offset rows and over @p offset
   * columns. If this object is initialized with an IndexSet, local_lines are
   * shifted as well.
   *
   * This function is useful if you are building block matrices, where all
   * blocks are built by the same DoFHandler object, i.e. the matrix size is
   * larger than the number of degrees of freedom. Since several matrix rows
   * and columns correspond to the same degrees of freedom, you'd generate
   * several constraint objects, then shift them, and finally merge() them
   * together again.
   */
  void
  shift(const size_type offset);

  /**
   * Clear all entries of this matrix. Reset the flag determining whether new
   * entries are accepted or not.
   *
   * This function may be called also on objects which are empty or already
   * cleared.
   */
  void
  clear();

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
  size_type
  n_constraints() const;

  /**
   * Return number of constraints stored in this matrix that are identities,
   * i.e., those constraints with only one degree of freedom and weight one.
   */
  size_type
  n_identities() const;

  /**
   * Return number of constraints stored in this matrix that have an
   * inhomogenity, i.e., those constraints with a non-trivial inhomogeneous
   * value.
   */
  size_type
  n_inhomogeneities() const;

  /**
   * Return whether the degree of freedom with number @p line_n is a
   * constrained one.
   *
   * Note that if close() was called before, then this function is
   * significantly faster, since then the constrained degrees of freedom are
   * sorted and we can do a binary search, while before close() was called, we
   * have to perform a linear search through all entries.
   */
  bool
  is_constrained(const size_type line_n) const;

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
  bool
  is_identity_constrained(const size_type line_n) const;

  /**
   * Return whether the two given degrees of freedom are linked by an equality
   * constraint that either constrains index1 to be so that
   * <code>index1=index2</code> or constrains index2 so that
   * <code>index2=index1</code>.
   */
  bool
  are_identity_constrained(const size_type line_n_1,
                           const size_type line_n_2) const;

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
  size_type
  max_constraint_indirections() const;

  /**
   * Return <tt>true</tt> in case the dof is constrained and there is a non-
   * trivial inhomogeneous values set to the dof.
   */
  bool
  is_inhomogeneously_constrained(const size_type index) const;

  /**
   * Return <tt>false</tt> if all constraints in the AffineConstraints are
   * homogeneous ones, and <tt>true</tt> if there is at least one
   * inhomogeneity.
   */
  bool
  has_inhomogeneities() const;

  /**
   * Return a pointer to the vector of entries if a line is constrained,
   * and a zero pointer in case the dof is not constrained.
   */
  const std::vector<std::pair<size_type, number>> *
  get_constraint_entries(const size_type line_n) const;

  /**
   * Return the value of the inhomogeneity stored in the constrained dof @p
   * line_n. Unconstrained dofs also return a zero value.
   */
  number
  get_inhomogeneity(const size_type line_n) const;

  /**
   * Print the constraints represented by the current object to the
   * given stream.
   *
   * For each constraint of the form
   * @f[
   *  x_{42} = 0.5 x_2 + 0.25 x_{14} + 2.75
   * @f]
   * this function will write a sequence of lines that look like this:
   * @code
   *   42 2 : 0.5
   *   42 14 : 0.25
   *   42 : 2.75
   * @endcode
   * The last line is only shown if the inhomogeneity (here: 2.75) is
   * nonzero.
   *
   * A block of lines such as the one above is repeated for each
   * constrained degree of freedom.
   */
  void
  print(std::ostream &out) const;

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
  void
  write_dot(std::ostream &) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Add the constraint indices associated to the indices in the given vector.
   * After a call to this function, the indices vector contains the initial
   * elements and all the associated constrained indices. This function sorts
   * the elements and suppresses duplicates.
   */
  void
  resolve_indices(std::vector<types::global_dof_index> &indices) const;

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
   * given sparsity pattern must not be compressed. The current object must be
   * closed. The sparsity pattern is compressed at the end of the function.
   */
  void
  condense(SparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square block sparsity patterns.
   */
  void
  condense(BlockSparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square compressed sparsity
   * patterns.
   */
  void
  condense(DynamicSparsityPattern &sparsity) const;

  /**
   * Same function as above, but condenses square compressed sparsity
   * patterns.
   */
  void
  condense(BlockDynamicSparsityPattern &sparsity) const;

  /**
   * Condense a given matrix, i.e., eliminate the rows and columns of the
   * matrix that correspond to constrained degrees of freedom.
   *
   * See the general documentation of this class for more detailed
   * information.
   */
  void
  condense(SparseMatrix<number> &matrix) const;

  /**
   * Same function as above, but condenses square block sparse matrices.
   */
  void
  condense(BlockSparseMatrix<number> &matrix) const;

  /**
   * Condense the given vector in-place. The @p VectorType may be a
   * Vector<float>, Vector<number>, BlockVector<tt><...></tt>, a PETSc or
   * Trilinos vector wrapper class, or any other type having the same
   * interface. Note that this function does not take any inhomogeneity into
   * account and throws an exception in case there are any inhomogeneities.
   * Use the function using both a matrix and vector for that case.
   *
   * @note This function does not work for MPI vectors. Use condense() with
   * two vector arguments instead.
   */
  template <typename VectorType>
  void
  condense(VectorType &vec) const;

  /**
   * The function copies and condenses values from @p vec_ghosted into @p
   * output. In a serial code it is equivalent to calling condense (vec). If
   * called in parallel, @p vec_ghosted is supposed to contain ghost elements
   * while @p output should not.
   */
  template <typename VectorType>
  void
  condense(const VectorType &vec_ghosted, VectorType &output) const;

  /**
   * Condense a given matrix and a given vector by eliminating rows and
   * columns of the linear system that correspond to constrained degrees of
   * freedom. The sparsity pattern associated with the matrix needs to be
   * condensed and compressed.  This function is the appropriate choice for
   * applying inhomogeneous constraints.
   *
   * The current object must be closed to call this function.
   *
   * See the general documentation of this class for more detailed
   * information.
   */
  template <typename VectorType>
  void
  condense(SparseMatrix<number> &matrix, VectorType &vector) const;

  /**
   * Same function as above, but condenses square block sparse matrices and
   * vectors.
   */
  template <typename BlockVectorType>
  void
  condense(BlockSparseMatrix<number> &matrix, BlockVectorType &vector) const;

  /**
   * Set the values of all constrained DoFs in a vector to zero.  The @p
   * VectorType may be a Vector<float>, Vector<number>,
   * BlockVector<tt><...></tt>, a PETSc or Trilinos vector wrapper class, or
   * any other type having the same interface.
   */
  template <typename VectorType>
  void
  set_zero(VectorType &vec) const;

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
   * local_dof_indices and distributes them to the global vector. In other
   * words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   * In most
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
  distribute_local_to_global(const InVector &              local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             OutVector &                   global_vector) const;

  /**
   * This function takes a vector of local contributions (@p local_vector)
   * corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global vector. In other
   * words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   * In most
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
   * stepping algorithm where the @ref GlossStiffnessMatrix "stiffness matrix" is assembled once, and the
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
  template <typename VectorType>
  void
  distribute_local_to_global(const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             VectorType &                  global_vector,
                             const FullMatrix<number> &    local_matrix) const;

  /**
   * Same as the previous function, except that it uses two (possibly) different
   * index sets to correctly handle inhomogeneities when the local matrix is
   * computed from a combination of two neighboring elements, for example for an
   * edge integral term in DG. Note that in the case that these two elements
   * have different polynomial degree, the local matrix is rectangular.
   *
   * <tt>local_dof_indices_row</tt> is the set of row indices and
   * <tt>local_dof_indices_col</tt> is the set of column indices of the local
   * matrix. <tt>diagonal=false</tt> says whether the two index sets are equal
   * or not.
   *
   * If both index sets are equal, <tt>diagonal</tt> must be set to true or we
   * simply use the previous function. If both index sets are different
   * (diagonal=false) the <tt>global_vector</tt> is modified to handle
   * inhomogeneities but no entries from <tt>local_vector</tt> are added. Note
   * that the edge integrals for inner edged for DG do not contribute any values
   * to the right hand side.
   */
  template <typename VectorType>
  void
  distribute_local_to_global(
    const Vector<number> &        local_vector,
    const std::vector<size_type> &local_dof_indices_row,
    const std::vector<size_type> &local_dof_indices_col,
    VectorType &                  global_vector,
    const FullMatrix<number> &    local_matrix,
    bool                          diagonal = false) const;

  /**
   * Enter a single value into a result vector, obeying constraints.
   */
  template <typename VectorType>
  void
  distribute_local_to_global(const size_type index,
                             const number    value,
                             VectorType &    global_vector) const;

  /**
   * This function takes a pointer to a vector of local contributions (@p
   * local_vector) corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global vector. In other
   * words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   * In most
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
  template <typename ForwardIteratorVec,
            typename ForwardIteratorInd,
            typename VectorType>
  void
  distribute_local_to_global(ForwardIteratorVec local_vector_begin,
                             ForwardIteratorVec local_vector_end,
                             ForwardIteratorInd local_indices_begin,
                             VectorType &       global_vector) const;

  /**
   * This function takes a matrix of local contributions (@p local_matrix)
   * corresponding to the degrees of freedom indices given in @p
   * local_dof_indices and distributes them to the global matrix. In other
   * words, this function implements a
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)).
   * In most
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
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * This function does almost the same as the function above but can treat
   * general rectangular matrices. The main difference to achieve this is that
   * the diagonal entries in constrained rows are left untouched instead of
   * being filled with arbitrary values.
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
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &row_indices,
                             const std::vector<size_type> &col_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * This function does almost the same as the function above for general
   * rectangular matrices but uses different AffineConstraints objects on the
   * row and column indices. The convention is that row indices are constrained
   * according to the calling AffineConstraints <code>*this</code>, whereas
   * column indices are constrained according to the given AffineConstraints
   * <code>column_affine_constraints</code>. This function allows to handle the
   * case where rows and columns of a matrix are represented by different
   * function spaces with their own enumeration of indices, as e.g. in mixed
   * finite element problems with separate DoFHandler objects or for flux
   * matrices between different levels in multigrid methods.
   *
   * Like the other method with separate slots for row and column indices,
   * this method does not add diagonal entries to eliminated degrees of
   * freedom. See there for a more elaborate description.
   */
  template <typename MatrixType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &row_indices,
                             const AffineConstraints &column_affine_constraints,
                             const std::vector<size_type> &column_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * This function simultaneously writes elements into matrix and vector,
   * according to the constraints specified by the calling AffineConstraints.
   * In other words, it performs the
   * [scatter
   * operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))
   * of the corresponding functions for matrices and vectors at the same time.
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
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             bool use_inhomogeneities_for_rhs = false) const;

  /**
   * Do a similar operation as the distribute_local_to_global() function that
   * distributes writing entries into a matrix for constrained degrees of
   * freedom, except that here we don't write into a matrix but only allocate
   * sparsity pattern entries.
   *
   * As explained in the
   * @ref hp_paper "hp-paper"
   * and in step-27, first allocating a sparsity pattern and later coming back
   * and allocating additional entries for those matrix entries that will be
   * written to due to the elimination of constrained degrees of freedom
   * (using AffineConstraints::condense() ), can be a very expensive procedure.
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
   * DoFTools::make_sparsity_pattern() function when passed an
   * AffineConstraints object.
   *
   * @note This function in itself is thread-safe, i.e., it works properly
   * also when several threads call it simultaneously. However, the function
   * call is only thread-safe if the underlying global sparsity pattern allows
   * for simultaneous access and the access is not to rows with the same
   * global index at the same time. This needs to be made sure from the
   * caller's site. There is no locking mechanism inside this method to
   * prevent data races.
   */
  void
  add_entries_local_to_global(
    const std::vector<size_type> &local_dof_indices,
    SparsityPatternBase &         sparsity_pattern,
    const bool                    keep_constrained_entries = true,
    const Table<2, bool> &        dof_mask = Table<2, bool>()) const;

  /**
   * Similar to the other function, but for non-quadratic sparsity patterns.
   */
  void
  add_entries_local_to_global(
    const std::vector<size_type> &row_indices,
    const std::vector<size_type> &col_indices,
    SparsityPatternBase &         sparsity_pattern,
    const bool                    keep_constrained_entries = true,
    const Table<2, bool> &        dof_mask = Table<2, bool>()) const;

  /**
   * Similar to the other function, but for non-quadratic sparsity patterns, and
   * for different constraints in the column space.
   */
  void
  add_entries_local_to_global(
    const std::vector<size_type> &   row_indices,
    const AffineConstraints<number> &col_constraints,
    const std::vector<size_type> &   col_indices,
    SparsityPatternBase &            sparsity_pattern,
    const bool                       keep_constrained_entries = true,
    const Table<2, bool> &           dof_mask = Table<2, bool>()) const;

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
  template <typename ForwardIteratorVec,
            typename ForwardIteratorInd,
            typename VectorType>
  void
  get_dof_values(const VectorType & global_vector,
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
   * Given a vector, set all constrained degrees of freedom to values so
   * that the constraints are satisfied. For example, if the current object
   * stores the constraint $x_3=\frac 12 x_1 + \frac 12 x_2$, then this
   * function will read the values of $x_1$ and $x_2$ from the given vector
   * and set the element $x_3$ according to this constraints. Similarly, if
   * the current object stores the constraint $x_{42}=208$, then this
   * function will set the 42nd element of the given vector to 208.
   *
   * @note If this function is called with a parallel vector @p vec, then the
   * vector must not contain ghost elements.
   */
  template <typename VectorType>
  void
  distribute(VectorType &vec) const;

  /**
   * @}
   */

  /**
   * This class represents one constraint in an AffineConstraints object.
   */
  struct ConstraintLine
  {
    /**
     * A data type in which we store the list of entries that make up the
     * homogeneous part of a constraint.
     */
    using Entries = std::vector<std::pair<size_type, number>>;

    /**
     * Global DoF index of this line. Since only very few lines are stored,
     * we can not assume a specific order and have to store the index
     * explicitly.
     */
    size_type index;

    /**
     * Row numbers and values of the entries in this line.
     *
     * For the reason why we use a vector instead of a map and the
     * consequences thereof, the same applies as what is said for
     * AffineConstraints::lines.
     */
    Entries entries;

    /**
     * Value of the inhomogeneity.
     */
    number inhomogeneity;

    /**
     * Default constructor.
     */
    ConstraintLine(const size_type &index = numbers::invalid_dof_index,
                   const typename AffineConstraints<
                     number>::ConstraintLine::Entries &entries       = {},
                   const number                        inhomogeneity = 0.0);

    /**
     * Copy constructor.
     */
    ConstraintLine(const ConstraintLine &other) = default;

    /**
     * Move constructor.
     */
    ConstraintLine(ConstraintLine &&other) noexcept = default;

    /**
     * Copy assignment.
     */
    ConstraintLine &
    operator=(const ConstraintLine &other) = default;

    /**
     * Move assignment.
     */
    ConstraintLine &
    operator=(ConstraintLine &&other) noexcept = default;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Write and read the data of this object from a stream for the purpose
     * of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int)
    {
      ar &index &entries &inhomogeneity;
    }

    /**
     * Swap function.
     */
    friend void
    swap(ConstraintLine &l1, ConstraintLine &l2)
    {
      std::swap(l1.index, l2.index);
      std::swap(l1.entries, l2.entries);
      std::swap(l1.inhomogeneity, l2.inhomogeneity);
    }
  };

  /**
   * Alias for the iterator type that is used in the LineRange container.
   */
  using const_iterator = typename std::vector<ConstraintLine>::const_iterator;

  /**
   * Alias for the return type used by get_lines().
   */
  using LineRange = boost::iterator_range<const_iterator>;

  /**
   * Return a range object containing (const) iterators to all line entries
   * stored in the AffineConstraints container. Such a range is useful to
   * initialize range-based for loops as supported by C++11.
   *
   * @return A range object for the half open range <code>[this->begin(),
   * this->end())</code> of line entries.
   */
  LineRange
  get_lines() const;

  /**
   * Check if the current object is consistent on all processors
   * in a distributed computation.
   *
   * This method checks if all processors agree on the constraints for their
   * local lines as given by @p locally_active_dofs. This method is a collective
   * operation and will return @p true only if all processors are consistent.
   *
   * Please supply the owned DoFs per processor as returned by
   * Utilities::MPI::all_gather(MPI_Comm, DoFHandler::locally_owned_dofs()) as
   * @p locally_owned_dofs and the result of
   * DoFTools::extract_locally_active_dofs() as
   * @p locally_active_dofs. The former is used to determine ownership of the
   * specific DoF, while the latter is used as the set of rows that need to be
   * checked.
   *
   * If @p verbose is set to @p true, additional debug information is written
   * to std::cout.
   *
   * @note This method exchanges all constraint information of locally active
   * lines and is as such slow for large computations and should probably
   * only be used in debug mode. We do not check all lines returned by
   * get_local_lines() but only the locally active ones, as we allow processors
   * to not know about some locally relevant rows.
   *
   * @return Whether all AffineConstraints objects are consistent. Returns
   * the same value on all processors.
   */
  bool
  is_consistent_in_parallel(const std::vector<IndexSet> &locally_owned_dofs,
                            const IndexSet &             locally_active_dofs,
                            const MPI_Comm               mpi_communicator,
                            const bool                   verbose = false) const;

  /**
   * Make the current object consistent on all processors
   * in a distributed computation. One should call this function before
   * calling close().
   */
  void
  make_consistent_in_parallel(const IndexSet &locally_owned_dofs,
                              const IndexSet &locally_relevant_dofs,
                              const MPI_Comm  mpi_communicator);

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcMatrixIsClosed);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcMatrixNotClosed);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1(ExcLineInexistant,
                 size_type,
                 << "The specified line " << arg1 << " does not exist.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException4(ExcEntryAlreadyExists,
                 size_type,
                 size_type,
                 number,
                 number,
                 << "The entry for the indices " << arg1 << " and " << arg2
                 << " already exists, but the values " << arg3 << " (old) and "
                 << arg4 << " (new) differ "
                 << "by " << (arg4 - arg3) << '.');
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2(ExcDoFConstrainedToConstrainedDoF,
                 int,
                 int,
                 << "You tried to constrain DoF " << arg1 << " to DoF " << arg2
                 << ", but that one is also constrained. This is not allowed!");
  /**
   * Exception.
   *
   * @ingroup Exceptions
   */
  DeclException1(ExcDoFIsConstrainedFromBothObjects,
                 size_type,
                 << "Degree of freedom " << arg1
                 << " is constrained from both object in a merge operation.");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1(ExcDoFIsConstrainedToConstrainedDoF,
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
  DeclException1(ExcRowNotStoredHere,
                 size_type,
                 << "The index set given to this constraints object indicates "
                 << "constraints for degree of freedom " << arg1
                 << " should not be stored by this object, but a constraint "
                 << "is being added.");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2(ExcColumnNotStoredHere,
                 size_type,
                 size_type,
                 << "The index set given to this constraints object indicates "
                 << "constraints using degree of freedom " << arg2
                 << " should not be stored by this object, but a constraint "
                 << "for degree of freedom " << arg1 << " uses it.");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException2(ExcIncorrectConstraint,
                 int,
                 int,
                 << "While distributing the constraint for DoF " << arg1
                 << ", it turns out that one of the processors "
                 << "who own the " << arg2 << " degrees of freedom that x_"
                 << arg1 << " is constrained against does not know about "
                 << "the constraint on x_" << arg1
                 << ". Did you not initialize the AffineConstraints container "
                 << "with the appropriate locally_relevant set so "
                 << "that every processor who owns a DoF that constrains "
                 << "another DoF also knows about this constraint?");

  template <typename>
  friend class AffineConstraints;

private:
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
   * This IndexSet is used to limit the lines to save in the AffineConstraints
   * to a subset. This is necessary, because the lines_cache vector would
   * become too big in a distributed calculation.
   */
  IndexSet local_lines;

  /**
   * Store whether the arrays are sorted.  If so, no new entries can be added.
   */
  bool sorted;

  mutable Threads::ThreadLocalStorage<
    internal::AffineConstraints::ScratchData<number>>
    scratch_data;

  /**
   * Internal function to calculate the index of line @p line_n in the vector
   * lines_cache using local_lines.
   */
  size_type
  calculate_line_index(const size_type line_n) const;

  /**
   * This function actually implements the local_to_global function for
   * standard (non-block) matrices.
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             const bool use_inhomogeneities_for_rhs,
                             const std::integral_constant<bool, false>) const;

  /**
   * This function actually implements the local_to_global function for block
   * matrices.
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             const bool use_inhomogeneities_for_rhs,
                             const std::integral_constant<bool, true>) const;

  /**
   * Internal helper function for distribute_local_to_global function.
   *
   * Creates a list of affected global rows for distribution, including the
   * local rows where the entries come from. The list is sorted according to
   * the global row indices.
   */
  void
  make_sorted_row_list(const std::vector<size_type> &local_dof_indices,
                       internal::AffineConstraints::GlobalRowsFromLocal<number>
                         &global_rows) const;

  /**
   * Internal helper function for add_entries_local_to_global function.
   *
   * Creates a list of affected rows for distribution without any additional
   * information, otherwise similar to the other make_sorted_row_list()
   * function.
   */
  void
  make_sorted_row_list(const std::vector<size_type> &local_dof_indices,
                       std::vector<size_type> &      active_dofs) const;

  /**
   * Internal helper function for distribute_local_to_global function.
   */
  template <typename MatrixScalar, typename VectorScalar>
  typename ProductType<VectorScalar, MatrixScalar>::type
  resolve_vector_entry(
    const size_type                                                 i,
    const internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows,
    const Vector<VectorScalar> &    local_vector,
    const std::vector<size_type> &  local_dof_indices,
    const FullMatrix<MatrixScalar> &local_matrix) const;
};

/* ---------------- template and inline functions ----------------- */

template <typename number>
inline AffineConstraints<number>::AffineConstraints(
  const IndexSet &local_constraints)
  : lines()
  , local_lines(local_constraints)
  , sorted(false)
{
  // make sure the IndexSet is compressed. Otherwise this can lead to crashes
  // that are hard to find (only happen in release mode).
  // see tests/mpi/affine_constraints_crash_01
  local_lines.compress();
}

template <typename number>
inline AffineConstraints<number>::AffineConstraints(
  const AffineConstraints &affine_constraints)
  : Subscriptor()
  , lines(affine_constraints.lines)
  , lines_cache(affine_constraints.lines_cache)
  , local_lines(affine_constraints.local_lines)
  , sorted(affine_constraints.sorted)
{}

template <typename number>
inline void
AffineConstraints<number>::add_line(const size_type line_n)
{
  Assert(sorted == false, ExcMatrixIsClosed());

  // the following can happen when we compute with distributed meshes and dof
  // handlers and we constrain a degree of freedom whose number we don't have
  // locally. if we don't abort here the program will try to allocate several
  // terabytes of memory to resize the various arrays below :-)
  Assert(line_n != numbers::invalid_size_type, ExcInternalError());
  const size_type line_index = calculate_line_index(line_n);

  // check whether line already exists; it may, in which case we can just quit
  if (is_constrained(line_n))
    return;

  // if necessary enlarge vector of existing entries for cache
  if (line_index >= lines_cache.size())
    lines_cache.resize(std::max(2 * static_cast<size_type>(lines_cache.size()),
                                line_index + 1),
                       numbers::invalid_size_type);

  // push a new line to the end of the list
  lines.emplace_back();
  lines.back().index         = line_n;
  lines.back().inhomogeneity = 0.;
  lines_cache[line_index]    = lines.size() - 1;
}



template <typename number>
inline void
AffineConstraints<number>::add_entry(const size_type constrained_dof_index,
                                     const size_type column,
                                     const number    weight)
{
  Assert(sorted == false, ExcMatrixIsClosed());
  Assert(constrained_dof_index != column,
         ExcMessage("Can't constrain a degree of freedom to itself"));

  // Ensure that the current line is present in the cache:
  const size_type line_index = calculate_line_index(constrained_dof_index);
  Assert(line_index < lines_cache.size(),
         ExcMessage("The current AffineConstraints does not contain the line "
                    "for the current entry. Call AffineConstraints::add_line "
                    "before calling this function."));

  // if in debug mode, check whether an entry for this column already exists
  // and if it's the same as the one entered at present
  //
  // in any case: exit the function if an entry for this column already
  // exists, since we don't want to enter it twice
  Assert(lines_cache[line_index] != numbers::invalid_size_type,
         ExcInternalError());
  Assert(!local_lines.size() || local_lines.is_element(column),
         ExcColumnNotStoredHere(constrained_dof_index, column));
  ConstraintLine *line_ptr = &lines[lines_cache[line_index]];
  Assert(line_ptr->index == constrained_dof_index, ExcInternalError());
  for (const auto &p : line_ptr->entries)
    if (p.first == column)
      {
        Assert(std::abs(p.second - weight) < 1.e-14,
               ExcEntryAlreadyExists(
                 constrained_dof_index, column, p.second, weight));
        return;
      }

  line_ptr->entries.emplace_back(column, weight);
}



template <typename number>
inline void
AffineConstraints<number>::set_inhomogeneity(
  const size_type constrained_dof_index,
  const number    value)
{
  const size_type line_index = calculate_line_index(constrained_dof_index);
  Assert(line_index < lines_cache.size() &&
           lines_cache[line_index] != numbers::invalid_size_type,
         ExcMessage("call add_line() before calling set_inhomogeneity()"));
  Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
  ConstraintLine *line_ptr = &lines[lines_cache[line_index]];
  line_ptr->inhomogeneity  = value;
}



template <typename number>
template <typename VectorType>
inline void
AffineConstraints<number>::set_zero(VectorType &vec) const
{
  // since lines is a private member, we cannot pass it to the functions
  // above. therefore, copy the content which is cheap
  std::vector<size_type> constrained_lines(lines.size());
  for (unsigned int i = 0; i < lines.size(); ++i)
    constrained_lines[i] = lines[i].index;
  internal::AffineConstraintsImplementation::set_zero_all(constrained_lines,
                                                          vec);
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::n_constraints() const
{
  return lines.size();
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::n_identities() const
{
  return std::count_if(lines.begin(),
                       lines.end(),
                       [](const ConstraintLine &line) {
                         return (line.entries.size() == 1) &&
                                (line.entries[0].second == number(1.));
                       });
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::n_inhomogeneities() const
{
  return std::count_if(lines.begin(),
                       lines.end(),
                       [](const ConstraintLine &line) {
                         return (line.inhomogeneity != number(0.));
                       });
}

template <typename number>
inline bool
AffineConstraints<number>::is_constrained(const size_type index) const
{
  if (lines.empty())
    return false;

  const size_type line_index = calculate_line_index(index);
  return ((line_index < lines_cache.size()) &&
          (lines_cache[line_index] != numbers::invalid_size_type));
}

template <typename number>
inline bool
AffineConstraints<number>::is_inhomogeneously_constrained(
  const size_type line_n) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return false;
  else
    {
      Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
      return (lines[lines_cache[line_index]].inhomogeneity != number(0.));
    }
}

template <typename number>
inline const std::vector<std::pair<types::global_dof_index, number>> *
AffineConstraints<number>::get_constraint_entries(const size_type line_n) const
{
  if (lines.empty())
    return nullptr;

  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return nullptr;
  else
    return &lines[lines_cache[line_index]].entries;
}

template <typename number>
inline number
AffineConstraints<number>::get_inhomogeneity(const size_type line_n) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return 0;
  else
    return lines[lines_cache[line_index]].inhomogeneity;
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::calculate_line_index(const size_type line_n) const
{
  // IndexSet is unused (serial case)
  if (local_lines.size() == 0)
    return line_n;

  Assert(local_lines.is_element(line_n), ExcRowNotStoredHere(line_n));

  return local_lines.index_within_set(line_n);
}

template <typename number>
inline bool
AffineConstraints<number>::can_store_line(size_type line_n) const
{
  return local_lines.size() == 0 || local_lines.is_element(line_n);
}

template <typename number>
inline const IndexSet &
AffineConstraints<number>::get_local_lines() const
{
  return local_lines;
}

template <typename number>
template <typename VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const size_type index,
  const number    value,
  VectorType &    global_vector) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());

  if (is_constrained(index) == false)
    global_vector(index) += value;
  else
    {
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(index)]];
      for (size_type j = 0; j < position.entries.size(); ++j)
        global_vector(position.entries[j].first) +=
          value * position.entries[j].second;
    }
}

template <typename number>
template <typename ForwardIteratorVec,
          typename ForwardIteratorInd,
          typename VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  ForwardIteratorVec local_vector_begin,
  ForwardIteratorVec local_vector_end,
  ForwardIteratorInd local_indices_begin,
  VectorType &       global_vector) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());
  for (; local_vector_begin != local_vector_end;
       ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        internal::ElementAccess<VectorType>::add(*local_vector_begin,
                                                 *local_indices_begin,
                                                 global_vector);
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          for (size_type j = 0; j < position.entries.size(); ++j)
            internal::ElementAccess<VectorType>::add(
              (*local_vector_begin) * position.entries[j].second,
              position.entries[j].first,
              global_vector);
        }
    }
}

template <typename number>
template <class InVector, class OutVector>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const InVector &              local_vector,
  const std::vector<size_type> &local_dof_indices,
  OutVector &                   global_vector) const
{
  Assert(local_vector.size() == local_dof_indices.size(),
         ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  distribute_local_to_global(local_vector.begin(),
                             local_vector.end(),
                             local_dof_indices.begin(),
                             global_vector);
}

template <typename number>
template <typename ForwardIteratorVec,
          typename ForwardIteratorInd,
          typename VectorType>
inline void
AffineConstraints<number>::get_dof_values(
  const VectorType & global_vector,
  ForwardIteratorInd local_indices_begin,
  ForwardIteratorVec local_vector_begin,
  ForwardIteratorVec local_vector_end) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());
  for (; local_vector_begin != local_vector_end;
       ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        *local_vector_begin = global_vector(*local_indices_begin);
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          typename VectorType::value_type value = position.inhomogeneity;
          for (size_type j = 0; j < position.entries.size(); ++j)
            value += (global_vector(position.entries[j].first) *
                      position.entries[j].second);
          *local_vector_begin = value;
        }
    }
}

// Forward declarations
#ifndef DOXYGEN
template <typename MatrixType>
class BlockMatrixBase;
template <typename SparsityPatternType>
class BlockSparsityPatternBase;
template <typename number>
class BlockSparseMatrixEZ;

namespace internal
{
  namespace AffineConstraints
  {
    /**
     * A "traits" class that can be used to determine whether a given type is a
     * block matrix type or not. For example,
     * @code
     *   IsBlockMatrix<SparseMatrix<number> >::value
     * @endcode
     * has the value `false`, whereas
     * @code
     *   IsBlockMatrix<BlockSparseMatrix<number> >::value
     * @endcode
     * is true. This is sometimes useful in template contexts where we may want
     * to do things differently depending on whether a template type denotes a
     * regular or a block matrix type.
     *
     * @see
     * @ref GlossBlockLA "Block (linear algebra)"
     */
    template <typename MatrixType>
    struct IsBlockMatrix
    {
    private:
      /**
       * Overload returning true if the class is derived from BlockMatrixBase,
       * which is what block matrices do (with the exception of
       * BlockSparseMatrixEZ).
       */
      template <typename T>
      static std::true_type
      check(const BlockMatrixBase<T> *);

      /**
       * Overload for BlockSparseMatrixEZ, which is the only block matrix not
       * derived from BlockMatrixBase at the time of writing this class.
       */
      template <typename T>
      static std::true_type
      check(const BlockSparseMatrixEZ<T> *);

      /**
       * Catch all for all other potential types that are then apparently not
       * block matrices.
       */
      static std::false_type
      check(...);

    public:
      /**
       * A statically computable value that indicates whether the template
       * argument to this class is a block matrix (in fact whether the type is
       * derived from BlockMatrixBase<T> or is one of the other block matrix
       * types).
       */
      static const bool value =
        std::is_same_v<decltype(check(std::declval<MatrixType *>())),
                       std::true_type>;
    };

    // instantiation of the static member
    template <typename MatrixType>
    const bool IsBlockMatrix<MatrixType>::value;

  } // namespace AffineConstraints
} // namespace internal
#endif



template <typename number>
template <typename other_number>
inline void
AffineConstraints<number>::copy_from(
  const AffineConstraints<other_number> &other)
{
  lines.clear();
  lines.reserve(other.lines.size());

  for (const auto l : other.lines)
    lines.emplace_back(l.index,
                       typename ConstraintLine::Entries(l.entries.begin(),
                                                        l.entries.end()),
                       l.inhomogeneity);

  lines_cache = other.lines_cache;
  local_lines = other.local_lines;
  sorted      = other.sorted;
}


template <typename number>
template <typename other_number>
void
AffineConstraints<number>::merge(
  const AffineConstraints<other_number> &other_constraints,
  const MergeConflictBehavior            merge_conflict_behavior,
  const bool                             allow_different_local_lines)
{
  (void)allow_different_local_lines;
  Assert(allow_different_local_lines ||
           local_lines == other_constraints.local_lines,
         ExcMessage(
           "local_lines for this and the other objects are not the same "
           "although allow_different_local_lines is false."));

  // store the previous state with respect to sorting
  const bool object_was_sorted = sorted;
  sorted                       = false;

  // first action is to fold into the present object possible constraints
  // in the second object. we don't strictly need to do this any more since
  // the AffineConstraints container has learned to deal with chains of
  // constraints in the close() function, but we have traditionally done
  // this and it's not overly hard to do.
  //
  // for this, loop over all constraints and replace the constraint lines
  // with a new one where constraints are replaced if necessary.
  typename ConstraintLine::Entries tmp;
  for (ConstraintLine &line : lines)
    {
      tmp.clear();
      for (const std::pair<size_type, number> &entry : line.entries)
        {
          // if the present dof is not stored, or not constrained, or if we
          // won't take the constraint from the other object, then simply copy
          // it over
          if ((other_constraints.local_lines.size() != 0. &&
               other_constraints.local_lines.is_element(entry.first) ==
                 false) ||
              other_constraints.is_constrained(entry.first) == false ||
              ((merge_conflict_behavior != right_object_wins) &&
               other_constraints.is_constrained(entry.first) &&
               this->is_constrained(entry.first)))
            tmp.push_back(entry);
          else
            // otherwise resolve further constraints by replacing the old
            // entry by a sequence of new entries taken from the other
            // object, but with multiplied weights
            {
              const auto *other_entries =
                other_constraints.get_constraint_entries(entry.first);
              Assert(other_entries != nullptr, ExcInternalError());

              const number weight = entry.second;

              for (const auto &other_entry : *other_entries)
                tmp.emplace_back(other_entry.first,
                                 other_entry.second * weight);

              line.inhomogeneity +=
                other_constraints.get_inhomogeneity(entry.first) * weight;
            }
        }
      // finally exchange old and newly resolved line
      line.entries.swap(tmp);
    }

  if (local_lines.size() != 0)
    local_lines.add_indices(other_constraints.local_lines);

  {
    // do not bother to resize the lines cache exactly since it is pretty
    // cheap to adjust it along the way.
    std::fill(lines_cache.begin(),
              lines_cache.end(),
              numbers::invalid_size_type);

    // reset lines_cache for our own constraints
    size_type index = 0;
    for (const ConstraintLine &line : lines)
      {
        const size_type local_line_no = calculate_line_index(line.index);
        if (local_line_no >= lines_cache.size())
          lines_cache.resize(local_line_no + 1, numbers::invalid_size_type);
        lines_cache[local_line_no] = index++;
      }

    // Add other_constraints to lines cache and our list of constraints
    for (const auto &line : other_constraints.lines)
      {
        const size_type local_line_no = calculate_line_index(line.index);
        if (local_line_no >= lines_cache.size())
          {
            lines_cache.resize(local_line_no + 1, numbers::invalid_size_type);
            lines.emplace_back(line.index,
                               typename ConstraintLine::Entries(
                                 line.entries.begin(), line.entries.end()),
                               line.inhomogeneity);
            lines_cache[local_line_no] = index++;
          }
        else if (lines_cache[local_line_no] == numbers::invalid_size_type)
          {
            // there are no constraints for that line yet
            lines.emplace_back(line.index,
                               typename ConstraintLine::Entries(
                                 line.entries.begin(), line.entries.end()),
                               line.inhomogeneity);
            AssertIndexRange(local_line_no, lines_cache.size());
            lines_cache[local_line_no] = index++;
          }
        else
          {
            // we already store that line
            switch (merge_conflict_behavior)
              {
                case no_conflicts_allowed:
                  AssertThrow(false,
                              ExcDoFIsConstrainedFromBothObjects(line.index));
                  break;

                case left_object_wins:
                  // ignore this constraint
                  break;

                case right_object_wins:
                  AssertIndexRange(local_line_no, lines_cache.size());
                  lines[lines_cache[local_line_no]] = {
                    line.index,
                    typename ConstraintLine::Entries(line.entries.begin(),
                                                     line.entries.end()),
                    static_cast<number>(line.inhomogeneity)};
                  break;

                default:
                  Assert(false, ExcNotImplemented());
              }
          }
      }

    // check that we set the pointers correctly
    for (size_type i = 0; i < lines_cache.size(); ++i)
      if (lines_cache[i] != numbers::invalid_size_type)
        Assert(i == calculate_line_index(lines[lines_cache[i]].index),
               ExcInternalError());
  }

  // if the object was sorted before, then make sure it is so afterward as
  // well. otherwise leave everything in the unsorted state
  if (object_was_sorted == true)
    close();
}



template <typename number>
template <typename MatrixType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number> &    local_matrix,
  const std::vector<size_type> &local_dof_indices,
  MatrixType &                  global_matrix) const
{
  // create a dummy and hand on to the function actually implementing this
  // feature in the cm.templates.h file.
  Vector<typename MatrixType::value_type> dummy(0);
  distribute_local_to_global(
    local_matrix,
    dummy,
    local_dof_indices,
    global_matrix,
    dummy,
    false,
    std::integral_constant<
      bool,
      internal::AffineConstraints::IsBlockMatrix<MatrixType>::value>());
}



template <typename number>
template <typename MatrixType, typename VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number> &    local_matrix,
  const Vector<number> &        local_vector,
  const std::vector<size_type> &local_dof_indices,
  MatrixType &                  global_matrix,
  VectorType &                  global_vector,
  bool                          use_inhomogeneities_for_rhs) const
{
  // enter the internal function with the respective block information set,
  // the actual implementation follows in the cm.templates.h file.
  distribute_local_to_global(
    local_matrix,
    local_vector,
    local_dof_indices,
    global_matrix,
    global_vector,
    use_inhomogeneities_for_rhs,
    std::integral_constant<
      bool,
      internal::AffineConstraints::IsBlockMatrix<MatrixType>::value>());
}



template <typename number>
inline AffineConstraints<number>::ConstraintLine::ConstraintLine(
  const size_type &                                                  index,
  const typename AffineConstraints<number>::ConstraintLine::Entries &entries,
  const number inhomogeneity)
  : index(index)
  , entries(entries)
  , inhomogeneity(inhomogeneity)
{}



DEAL_II_NAMESPACE_CLOSE

#endif
