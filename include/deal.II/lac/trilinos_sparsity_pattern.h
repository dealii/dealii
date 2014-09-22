// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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

#ifndef __deal2__trilinos_sparsity_pattern_h
#define __deal2__trilinos_sparsity_pattern_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/lac/exceptions.h>

#  include <vector>
#  include <cmath>
#  include <memory>

#  include <deal.II/base/std_cxx11/shared_ptr.h>

#  include <Epetra_FECrsGraph.h>
#  include <Epetra_Map.h>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#    include "mpi.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif

DEAL_II_NAMESPACE_OPEN

// forward declarations
class SparsityPattern;
class CompressedSparsityPattern;
class CompressedSetSparsityPattern;
class CompressedSimpleSparsityPattern;

namespace TrilinosWrappers
{
  // forward declarations
  class SparsityPattern;

  namespace SparsityPatternIterators
  {
    // forward declaration
    class Iterator;

    /**
     * Accessor class for iterators into sparsity patterns. This class is
     * also the base class for both const and non-const accessor classes
     * into sparse matrices.
     *
     * Note that this class only allows read access to elements, providing
     * their row and column number. It does not allow modifying the
     * sparsity pattern itself.
     *
     * @ingroup TrilinosWrappers
     * @author Wolfgang Bangerth, Martin Kronbichler, Guido Kanschat
     * @date 2004, 2008, 2012
     */
    class Accessor
    {
    public:
      /**
       * Declare type for container size.
       */
      typedef dealii::types::global_dof_index size_type;

      /**
       * Constructor.
       */
      Accessor (const SparsityPattern *sparsity_pattern,
                const size_type        row,
                const size_type        index);

      /**
       * Copy constructor.
       */
      Accessor (const Accessor &a);

      /**
       * Row number of the element represented by this object.
       */
      size_type row() const;

      /**
       * Index in row of the element represented by this object.
       */
      size_type index() const;

      /**
       * Column number of the element represented by this object.
       */
      size_type column() const;

      /**
       * Exception
       */
      DeclException0 (ExcBeyondEndOfSparsityPattern);

      /**
       * Exception
       */
      DeclException3 (ExcAccessToNonlocalRow,
                      size_type, size_type, size_type,
                      << "You tried to access row " << arg1
                      << " of a distributed sparsity pattern, "
                      << " but only rows " << arg2 << " through " << arg3
                      << " are stored locally and can be accessed.");

    private:
      /**
       * The matrix accessed.
       */
      mutable SparsityPattern *sparsity_pattern;

      /**
       * Current row number.
       */
      size_type a_row;

      /**
       * Current index in row.
       */
      size_type a_index;

      /**
       * Cache where we store the column indices of the present row. This is
       * necessary, since Trilinos makes access to the elements of its
       * matrices rather hard, and it is much more efficient to copy all
       * column entries of a row once when we enter it than repeatedly asking
       * Trilinos for individual ones. This also makes some sense since it is
       * likely that we will access them sequentially anyway.
       *
       * In order to make copying of iterators/accessor of acceptable
       * performance, we keep a shared pointer to these entries so that more
       * than one accessor can access this data if necessary.
       */
      std_cxx11::shared_ptr<const std::vector<size_type> > colnum_cache;

      /**
       * Discard the old row caches (they may still be used by other
       * accessors) and generate new ones for the row pointed to presently by
       * this accessor.
       */
      void visit_present_row ();

      /**
       * Make enclosing class a friend.
       */
      friend class Iterator;
    };

    /**
     * Iterator class for sparsity patterns of type TrilinosWrappers::SparsityPattern.
     * Access to individual elements of the sparsity pattern is handled by the
     * Accessor class in this namespace.
     */
    class Iterator
    {
    public:
      /**
       * Declare type for container size.
       */
      typedef dealii::types::global_dof_index size_type;

      /**
       * Constructor. Create an iterator into the matrix @p matrix for the
       * given row and the index within it.
       */
      Iterator (const SparsityPattern *sparsity_pattern,
                const size_type        row,
                const size_type        index);

      /**
       * Copy constructor.
       */
      Iterator (const Iterator &i);

      /**
       * Prefix increment.
       */
      Iterator &operator++ ();

      /**
       * Postfix increment.
       */
      Iterator operator++ (int);

      /**
       * Dereferencing operator.
       */
      const Accessor &operator* () const;

      /**
       * Dereferencing operator.
       */
      const Accessor *operator-> () const;

      /**
       * Comparison. True, if both iterators point to the same matrix
       * position.
       */
      bool operator == (const Iterator &) const;

      /**
       * Inverse of <tt>==</tt>.
       */
      bool operator != (const Iterator &) const;

      /**
       * Comparison operator. Result is true if either the first row number is
       * smaller or if the row numbers are equal and the first index is
       * smaller.
       */
      bool operator < (const Iterator &) const;

      /**
       * Exception
       */
      DeclException2 (ExcInvalidIndexWithinRow,
                      size_type, size_type,
                      << "Attempt to access element " << arg2
                      << " of row " << arg1
                      << " which doesn't have that many elements.");

    private:
      /**
       * Store an object of the accessor class.
       */
      Accessor accessor;

      friend class TrilinosWrappers::SparsityPattern;
    };

  }


  /**
   * This class implements a wrapper class to use the Trilinos distributed
   * sparsity pattern class Epetra_FECrsGraph. This class is designed to be
   * used for construction of %parallel Trilinos matrices. The functionality of
   * this class is modeled after the existing sparsity pattern classes, with
   * the difference that this class can work fully in %parallel according to a
   * partitioning of the sparsity pattern rows.
   *
   * This class has many similarities to the compressed sparsity pattern
   * classes of deal.II (i.e., the classes CompressedSparsityPattern,
   * CompressedSetSparsityPattern, and CompressedSimpleSparsityPattern), since
   * it can dynamically add elements to the pattern without any memory being
   * previously reserved for it. However, it also has a method
   * SparsityPattern::compress(), that finalizes the pattern and enables its
   * use with Trilinos sparse matrices.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Sparsity
   * @author Martin Kronbichler, 2008
   */
  class SparsityPattern : public Subscriptor
  {
  public:

    /**
     * Declare type for container size.
     */
    typedef dealii::types::global_dof_index size_type;

    /**
     * Declare a typedef for the
     * iterator class.
     */
    typedef SparsityPatternIterators::Iterator const_iterator;

    /**
     * @name Basic constructors and initalization.
     */
//@{
    /**
     * Default constructor. Generates an empty (zero-size) sparsity pattern.
     */
    SparsityPattern ();

    /**
     * Generate a sparsity pattern that is completely stored locally, having
     * $m$ rows and $n$ columns. The resulting matrix will be completely
     * stored locally, too.
     *
     * It is possible to specify the number of columns entries per row using
     * the optional @p n_entries_per_row argument. However, this value does
     * not need to be accurate or even given at all, since one does usually
     * not have this kind of information before building the sparsity pattern
     * (the usual case when the function DoFTools::make_sparsity_pattern() is
     * called). The entries are allocated dynamically in a similar manner as
     * for the deal.II CompressedSparsityPattern classes. However, a good
     * estimate will reduce the setup time of the sparsity pattern.
     */
    SparsityPattern (const size_type  m,
                     const size_type  n,
                     const size_type  n_entries_per_row = 0);

    /**
     * Generate a sparsity pattern that is completely stored locally, having
     * $m$ rows and $n$ columns. The resulting matrix will be completely
     * stored locally, too.
     *
     * The vector <tt>n_entries_per_row</tt> specifies the number of entries
     * in each row (an information usually not available, though).
     */
    SparsityPattern (const size_type               m,
                     const size_type               n,
                     const std::vector<size_type> &n_entries_per_row);

    /**
     * Copy constructor. Sets the calling sparsity pattern to be the same as
     * the input sparsity pattern.
     */
    SparsityPattern (const SparsityPattern &input_sparsity_pattern);

    /**
     * Destructor. Made virtual so that one can use pointers to this class.
     */
    virtual ~SparsityPattern ();

    /**
     * Initialize a sparsity pattern that is completely stored locally, having
     * $m$ rows and $n$ columns. The resulting matrix will be completely
     * stored locally.
     *
     * The number of columns entries per row is specified as the maximum
     * number of entries argument.  This does not need to be an accurate
     * number since the entries are allocated dynamically in a similar manner
     * as for the deal.II CompressedSparsityPattern classes, but a good
     * estimate will reduce the setup time of the sparsity pattern.
     */
    void
    reinit (const size_type  m,
            const size_type  n,
            const size_type  n_entries_per_row = 0);

    /**
     * Initialize a sparsity pattern that is completely stored locally, having
     * $m$ rows and $n$ columns. The resulting matrix will be completely
     * stored locally.
     *
     * The vector <tt>n_entries_per_row</tt> specifies the number of entries
     * in each row.
     */
    void
    reinit (const size_type               m,
            const size_type               n,
            const std::vector<size_type> &n_entries_per_row);

    /**
     * Copy function. Sets the calling sparsity pattern to be the same as the
     * input sparsity pattern.
     */
    void
    copy_from (const SparsityPattern &input_sparsity_pattern);

    /**
     * Copy function from one of the deal.II sparsity patterns. If used in
     * parallel, this function uses an ad-hoc partitioning of the rows and
     * columns.
     */
    template<typename SparsityType>
    void
    copy_from (const SparsityType &nontrilinos_sparsity_pattern);

    /**
     * Copy operator. This operation is only allowed for empty objects, to
     * avoid potentially very costly operations automatically synthesized by
     * the compiler. Use copy_from() instead if you know that you really want
     * to copy a sparsity pattern with non-trivial content.
     */
    SparsityPattern &operator = (const SparsityPattern &input_sparsity_pattern);

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     *
     * This is a collective operation that needs to be called on all
     * processors in order to avoid a dead lock.
     */
    void clear ();

    /**
     * In analogy to our own SparsityPattern class, this function compresses
     * the sparsity pattern and allows the resulting pattern to be used for
     * actually generating a (Trilinos-based) matrix. This function also
     * exchanges non-local data that might have accumulated during the
     * addition of new elements. This function must therefore be called once
     * the structure is fixed. This is a collective operation, i.e., it needs
     * to be run on all processors when used in parallel.
     */
    void compress ();
//@}
    /**
     * @name Constructors and initialization using an Epetra_Map description
     */
//@{

    /**
     * Constructor for a square sparsity pattern using an Epetra_map for the
     * description of the %parallel partitioning. Moreover, the number of
     * nonzero entries in the rows of the sparsity pattern can be
     * specified. Note that this number does not need to be exact, and it is
     * allowed that the actual sparsity structure has more nonzero entries
     * than specified in the constructor (the usual case when the function
     * DoFTools::make_sparsity_pattern() is called). However it is still
     * advantageous to provide good estimates here since a good value will
     * avoid repeated allocation of memory, which considerably increases the
     * performance when creating the sparsity pattern.
     */
    SparsityPattern (const Epetra_Map &parallel_partitioning,
                     const size_type   n_entries_per_row = 0);

    /**
     * Same as before, but now use the exact number of nonzeros in each m
     * row. Since we know the number of elements in the sparsity pattern
     * exactly in this case, we can already allocate the right amount of
     * memory, which makes the creation process by the respective
     * SparsityPattern::reinit call considerably faster. However, this is a
     * rather unusual situation, since knowing the number of entries in each
     * row is usually connected to knowing the indices of nonzero entries,
     * which the sparsity pattern is designed to describe.
     */
    SparsityPattern (const Epetra_Map             &parallel_partitioning,
                     const std::vector<size_type> &n_entries_per_row);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular sparsity pattern, where one map
     * describes the %parallel partitioning of the dofs associated with the
     * sparsity pattern rows and the other one of the sparsity pattern
     * columns. Note that there is no real parallelism along the columns
     * &ndash; the processor that owns a certain row always owns all the
     * column elements, no matter how far they might be spread out. The second
     * Epetra_Map is only used to specify the number of columns and for
     * specifying the correct domain space when performing matrix-vector
     * products with vectors based on the same column map.
     *
     * The number of columns entries per row is specified as the maximum
     * number of entries argument.
     */
    SparsityPattern (const Epetra_Map   &row_parallel_partitioning,
                     const Epetra_Map   &col_parallel_partitioning,
                     const size_type     n_entries_per_row = 0);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular matrices, where one map specifies
     * the %parallel distribution of rows and the second one specifies the
     * distribution of degrees of freedom associated with matrix columns. This
     * second map is however not used for the distribution of the columns
     * themselves &ndash; rather, all column elements of a row are stored on
     * the same processor. The vector <tt>n_entries_per_row</tt> specifies the
     * number of entries in each row of the newly generated matrix.
     */
    SparsityPattern (const Epetra_Map             &row_parallel_partitioning,
                     const Epetra_Map             &col_parallel_partitioning,
                     const std::vector<size_type> &n_entries_per_row);

    /**
     * Reinitialization function for generating a square sparsity pattern
     * using an Epetra_Map for the description of the %parallel partitioning
     * and the number of nonzero entries in the rows of the sparsity
     * pattern. Note that this number does not need to be exact, and it is
     * even allowed that the actual sparsity structure has more nonzero
     * entries than specified in the constructor. However it is still
     * advantageous to provide good estimates here since this will
     * considerably increase the performance when creating the sparsity
     * pattern.
     *
     * This function does not create any entries by itself, but provides the
     * correct data structures that can be used by the respective add()
     * function.
     */
    void
    reinit (const Epetra_Map &parallel_partitioning,
            const size_type   n_entries_per_row = 0);

    /**
     * Same as before, but now use the exact number of nonzeros in each m
     * row. Since we know the number of elements in the sparsity pattern
     * exactly in this case, we can already allocate the right amount of
     * memory, which makes process of adding entries to the sparsity pattern
     * considerably faster. However, this is a rather unusual situation, since
     * knowing the number of entries in each row is usually connected to
     * knowing the indices of nonzero entries, which the sparsity pattern is
     * designed to describe.
     */
    void
    reinit (const Epetra_Map             &parallel_partitioning,
            const std::vector<size_type> &n_entries_per_row);

    /**
     * This reinit function is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular sparsity pattern, where one map
     * describes the %parallel partitioning of the dofs associated with the
     * sparsity pattern rows and the other one of the sparsity pattern
     * columns. Note that there is no real parallelism along the columns
     * &ndash; the processor that owns a certain row always owns all the
     * column elements, no matter how far they might be spread out. The second
     * Epetra_Map is only used to specify the number of columns and for
     * internal arragements when doing matrix-vector products with vectors
     * based on that column map.
     *
     * The number of columns entries per row is specified by the argument
     * <tt>n_entries_per_row</tt>.
     */
    void
    reinit (const Epetra_Map   &row_parallel_partitioning,
            const Epetra_Map   &col_parallel_partitioning,
            const size_type     n_entries_per_row = 0);

    /**
     * This reinit function is similar to the one above, but it now takes two
     * different Epetra maps for rows and columns. This interface is meant to
     * be used for generating rectangular matrices, where one map specifies
     * the %parallel distribution of rows and the second one specifies the
     * distribution of degrees of freedom associated with matrix columns. This
     * second map is however not used for the distribution of the columns
     * themselves &ndash; rather, all column elements of a row are stored on
     * the same processor. The vector <tt>n_entries_per_row</tt> specifies the
     * number of entries in each row of the newly generated matrix.
     */
    void
    reinit (const Epetra_Map             &row_parallel_partitioning,
            const Epetra_Map             &col_parallel_partitioning,
            const std::vector<size_type> &n_entries_per_row);

    /**
     * Reinit function. Takes one of the deal.II sparsity patterns and a
     * %parallel partitioning of the rows and columns for initializing the
     * current Trilinos sparsity pattern. The optional argument @p
     * exchange_data can be used for reinitialization with a sparsity pattern
     * that is not fully constructed. This feature is only implemented for
     * input sparsity patterns of type CompressedSimpleSparsityPattern.
     */
    template<typename SparsityType>
    void
    reinit (const Epetra_Map   &row_parallel_partitioning,
            const Epetra_Map   &col_parallel_partitioning,
            const SparsityType &nontrilinos_sparsity_pattern,
            const bool          exchange_data = false);

    /**
     * Reinit function. Takes one of the deal.II sparsity patterns and a
     * %parallel partitioning of the rows and columns for initializing the
     * current Trilinos sparsity pattern. The optional argument @p
     * exchange_data can be used for reinitialization with a sparsity pattern
     * that is not fully constructed. This feature is only implemented for
     * input sparsity patterns of type CompressedSimpleSparsityPattern.
     */
    template<typename SparsityType>
    void
    reinit (const Epetra_Map   &parallel_partitioning,
            const SparsityType &nontrilinos_sparsity_pattern,
            const bool          exchange_data = false);
//@}
    /**
     * @name Constructors and initialization using an IndexSet description
     */
//@{

    /**
     * Constructor for a square sparsity pattern using an IndexSet and an MPI
     * communicator for the description of the %parallel
     * partitioning. Moreover, the number of nonzero entries in the rows of
     * the sparsity pattern can be specified. Note that this number does not
     * need to be exact, and it is even allowed that the actual sparsity
     * structure has more nonzero entries than specified in the
     * constructor. However it is still advantageous to provide good estimates
     * here since a good value will avoid repeated allocation of memory, which
     * considerably increases the performance when creating the sparsity
     * pattern.
     */
    SparsityPattern (const IndexSet  &parallel_partitioning,
                     const MPI_Comm  &communicator = MPI_COMM_WORLD,
                     const size_type  n_entries_per_row = 0);

    /**
     * Same as before, but now use the exact number of nonzeros in each m
     * row. Since we know the number of elements in the sparsity pattern
     * exactly in this case, we can already allocate the right amount of
     * memory, which makes the creation process by the respective
     * SparsityPattern::reinit call considerably faster. However, this is a
     * rather unusual situation, since knowing the number of entries in each
     * row is usually connected to knowing the indices of nonzero entries,
     * which the sparsity pattern is designed to describe.
     */
    SparsityPattern (const IndexSet                  &parallel_partitioning,
                     const MPI_Comm                  &communicator,
                     const std::vector<size_type> &n_entries_per_row);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different index sets to describe the %parallel partitioning of rows and
     * columns. This interface is meant to be used for generating rectangular
     * sparsity pattern. Note that there is no real parallelism along the
     * columns &ndash; the processor that owns a certain row always owns all
     * the column elements, no matter how far they might be spread out. The
     * second Epetra_Map is only used to specify the number of columns and for
     * internal arragements when doing matrix-vector products with vectors
     * based on that column map.
     *
     * The number of columns entries per row is specified as the maximum
     * number of entries argument.
     */
    SparsityPattern (const IndexSet  &row_parallel_partitioning,
                     const IndexSet  &col_parallel_partitioning,
                     const MPI_Comm  &communicator = MPI_COMM_WORLD,
                     const size_type  n_entries_per_row = 0);

    /**
     * This constructor is similar to the one above, but it now takes two
     * different index sets for rows and columns. This interface is meant to
     * be used for generating rectangular matrices, where one map specifies
     * the %parallel distribution of rows and the second one specifies the
     * distribution of degrees of freedom associated with matrix columns. This
     * second map is however not used for the distribution of the columns
     * themselves &ndash; rather, all column elements of a row are stored on
     * the same processor. The vector <tt>n_entries_per_row</tt> specifies the
     * number of entries in each row of the newly generated matrix.
     */
    SparsityPattern (const IndexSet               &row_parallel_partitioning,
                     const IndexSet               &col_parallel_partitioning,
                     const MPI_Comm               &communicator,
                     const std::vector<size_type> &n_entries_per_row);

    /**
     * This constructor constructs general sparsity patterns, possible
     * non-square ones. Constructing a sparsity pattern this way allows the
     * user to explicitly specify the rows into which we are going to add
     * elements. This set is required to be a superset of the first index set
     * @p row_parallel_partitioning that includes also rows that are owned by
     * another processor (ghost rows). Note that elements can only be added to
     * rows specified by @p writable_rows.
     *
     * This method is beneficial when the rows to which a processor is going
     * to write can be determined before actually inserting elements into the
     * matrix. For the typical parallel::distributed::Triangulation class used
     * in deal.II, we know that a processor only will add row elements for
     * what we call the locally relevant dofs (see
     * DoFTools::extract_locally_relevant_dofs). The other constructors
     * methods use general Trilinos facilities that allow to add elements to
     * arbitrary rows (as done by all the other reinit functions). However,
     * this flexbility come at a cost, the most prominent being that adding
     * elements into the same matrix from multiple threads in shared memory is
     * not safe whenever MPI is used. For these settings, the current method
     * is the one to choose: It will store the off-processor data as an
     * additional sparsity pattern (that is then passed to the Trilinos matrix
     * via the reinit mehtod) which can be organized in such a way that
     * thread-safety can be ensured (as long as the user makes sure to never
     * write into the same matrix row simultaneously, of course).
     */
    SparsityPattern (const IndexSet  &row_parallel_partitioning,
                     const IndexSet  &col_parallel_partitioning,
                     const IndexSet  &writable_rows,
                     const MPI_Comm  &communicator = MPI_COMM_WORLD,
                     const size_type  n_entries_per_row = 0);

    /**
     * Reinitialization function for generating a square sparsity pattern
     * using an IndexSet and an MPI communicator for the description of the
     * %parallel partitioning and the number of nonzero entries in the rows of
     * the sparsity pattern. Note that this number does not need to be exact,
     * and it is even allowed that the actual sparsity structure has more
     * nonzero entries than specified in the constructor. However it is still
     * advantageous to provide good estimates here since this will
     * considerably increase the performance when creating the sparsity
     * pattern.
     *
     * This function does not create any entries by itself, but provides the
     * correct data structures that can be used by the respective add()
     * function.
     */
    void
    reinit (const IndexSet  &parallel_partitioning,
            const MPI_Comm  &communicator = MPI_COMM_WORLD,
            const size_type  n_entries_per_row = 0);

    /**
     * Same as before, but now use the exact number of nonzeros in each m
     * row. Since we know the number of elements in the sparsity pattern
     * exactly in this case, we can already allocate the right amount of
     * memory, which makes process of adding entries to the sparsity pattern
     * considerably faster. However, this is a rather unusual situation, since
     * knowing the number of entries in each row is usually connected to
     * knowing the indices of nonzero entries, which the sparsity pattern is
     * designed to describe.
     */
    void
    reinit (const IndexSet               &parallel_partitioning,
            const MPI_Comm               &communicator,
            const std::vector<size_type> &n_entries_per_row);

    /**
     * This reinit function is similar to the one above, but it now takes two
     * different index sets for rows and columns. This interface is meant to
     * be used for generating rectangular sparsity pattern, where one index
     * set describes the %parallel partitioning of the dofs associated with
     * the sparsity pattern rows and the other one of the sparsity pattern
     * columns. Note that there is no real parallelism along the columns
     * &ndash; the processor that owns a certain row always owns all the
     * column elements, no matter how far they might be spread out. The second
     * IndexSet is only used to specify the number of columns and for internal
     * arragements when doing matrix-vector products with vectors based on an
     * EpetraMap based on that IndexSet.
     *
     * The number of columns entries per row is specified by the argument
     * <tt>n_entries_per_row</tt>.
     */
    void
    reinit (const IndexSet  &row_parallel_partitioning,
            const IndexSet  &col_parallel_partitioning,
            const MPI_Comm  &communicator = MPI_COMM_WORLD,
            const size_type  n_entries_per_row = 0);

    /**
     * This reinit function is used to specify general matrices, possibly
     * non-square ones. In addition to the arguments of the other reinit
     * method above, it allows the user to explicitly specify the rows into
     * which we are going to add elements. This set is a superset of the first
     * index set @p row_parallel_partitioning that includes also rows that are
     * owned by another processor (ghost rows).
     *
     * This method is beneficial when the rows to which a processor is going
     * to write can be determined before actually inserting elements into the
     * matrix. For the typical parallel::distributed::Triangulation class used
     * in deal.II, we know that a processor only will add row elements for
     * what we call the locally relevant dofs (see
     * DoFTools::extract_locally_relevant_dofs). Trilinos matrices allow to
     * add elements to arbitrary rows (as done by all the other reinit
     * functions) and this is what all the other reinit methods do,
     * too. However, this flexbility come at a cost, the most prominent being
     * that adding elements into the same matrix from multiple threads in
     * shared memory is not safe whenever MPI is used. For these settings, the
     * current method is the one to choose: It will store the off-processor
     * data as an additional sparsity pattern (that is then passed to the
     * Trilinos matrix via the reinit mehtod) which can be organized in such a
     * way that thread-safety can be ensured (as long as the user makes sure
     * to never write into the same matrix row simultaneously, of course).
     */
    void
    reinit (const IndexSet  &row_parallel_partitioning,
            const IndexSet  &col_parallel_partitioning,
            const IndexSet  &writeable_rows,
            const MPI_Comm  &communicator = MPI_COMM_WORLD,
            const size_type  n_entries_per_row = 0);

    /**
     * Same as before, but now using a vector <tt>n_entries_per_row</tt> for
     * specifying the number of entries in each row of the sparsity pattern.
     */
    void
    reinit (const IndexSet               &row_parallel_partitioning,
            const IndexSet               &col_parallel_partitioning,
            const MPI_Comm               &communicator,
            const std::vector<size_type> &n_entries_per_row);

    /**
     * Reinit function. Takes one of the deal.II sparsity patterns and the
     * %parallel partitioning of the rows and columns specified by two index
     * sets and a %parallel communicator for initializing the current Trilinos
     * sparsity pattern. The optional argument @p exchange_data can be used
     * for reinitialization with a sparsity pattern that is not fully
     * constructed. This feature is only implemented for input sparsity
     * patterns of type CompressedSimpleSparsityPattern.
     */
    template<typename SparsityType>
    void
    reinit (const IndexSet     &row_parallel_partitioning,
            const IndexSet     &col_parallel_partitioning,
            const SparsityType &nontrilinos_sparsity_pattern,
            const MPI_Comm     &communicator = MPI_COMM_WORLD,
            const bool          exchange_data = false);

    /**
     * Reinit function. Takes one of the deal.II sparsity patterns and a
     * %parallel partitioning of the rows and columns for initializing the
     * current Trilinos sparsity pattern. The optional argument @p
     * exchange_data can be used for reinitialization with a sparsity pattern
     * that is not fully constructed. This feature is only implemented for
     * input sparsity patterns of type CompressedSimpleSparsityPattern.
     */
    template<typename SparsityType>
    void
    reinit (const IndexSet     &parallel_partitioning,
            const SparsityType &nontrilinos_sparsity_pattern,
            const MPI_Comm     &communicator = MPI_COMM_WORLD,
            const bool          exchange_data = false);
//@}
    /**
     * @name Information on the sparsity pattern
     */
//@{

    /**
     * Returns the state of the sparsity pattern, i.e., whether compress()
     * needs to be called after an operation requiring data exchange.
     */
    bool is_compressed () const;

    /**
     * Gives the maximum number of entries per row on the current processor.
     */
    unsigned int max_entries_per_row () const;

    /**
     * Return the number of rows in this sparsity pattern.
     */
    size_type n_rows () const;

    /**
     * Return the number of columns in this sparsity pattern.
     */
    size_type n_cols () const;

    /**
     * Return the local dimension of the sparsity pattern, i.e. the number of
     * rows stored on the present MPI process. In the sequential case, this
     * number is the same as n_rows(), but for parallel matrices it may be
     * smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range().
     */
    unsigned int local_size () const;

    /**
     * Return a pair of indices indicating which rows of this sparsity pattern
     * are stored locally. The first number is the index of the first row
     * stored, the second the index of the one past the last one that is
     * stored locally. If this is a sequential matrix, then the result will be
     * the pair (0,n_rows()), otherwise it will be a pair (i,i+n), where
     * <tt>n=local_size()</tt>.
     */
    std::pair<size_type, size_type>
    local_range () const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     */
    bool in_local_range (const size_type index) const;

    /**
     * Return the number of nonzero elements of this sparsity pattern.
     */
    size_type n_nonzero_elements () const;

    /**
     * Number of entries in a specific row.
     */
    size_type row_length (const size_type row) const;

    /**
     * Compute the bandwidth of the matrix represented by this structure. The
     * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
     * represents a nonzero entry of the matrix. Consequently, the maximum
     * bandwidth a $n\times m$ matrix can have is $\max\{n-1,m-1\}$.
     */
    size_type bandwidth () const;

    /**
     * Return whether the object is empty. It is empty if no memory is
     * allocated, which is the same as when both dimensions are zero.
     */
    bool empty () const;

    /**
     * Return whether the index (<i>i,j</i>) exists in the sparsity pattern
     * (i.e., it may be non-zero) or not.
     */
    bool exists (const size_type i,
                 const size_type j) const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object. Currently not implemented for this class.
     */
    std::size_t memory_consumption () const;

//@}
    /**
     * @name Adding entries
     */
//@{
    /**
     * Add the element (<i>i,j</i>) to the sparsity pattern.
     */
    void add (const size_type i,
              const size_type j);


    /**
     * Add several elements in one row to the sparsity pattern.
     */
    template <typename ForwardIterator>
    void add_entries (const size_type  row,
                      ForwardIterator  begin,
                      ForwardIterator  end,
                      const bool       indices_are_sorted = false);
//@}
    /**
     * @name Access of underlying Trilinos data
     */
//@{

    /**
     * Return a const reference to the underlying Trilinos Epetra_CrsGraph
     * data that stores the sparsity pattern.
     */
    const Epetra_FECrsGraph &trilinos_sparsity_pattern () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the parallel partitioning of the domain space of this sparsity
     * pattern, i.e., the partitioning of the vectors matrices based on this
     * sparsity pattern are multiplied with.
     */
    const Epetra_Map &domain_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the range space of this sparsity pattern,
     * i.e., the partitioning of the vectors that are result from
     * matrix-vector products.
     */
    const Epetra_Map &range_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the sparsity pattern rows. Equal to the
     * partitioning of the range.
     */
    const Epetra_Map &row_partitioner () const;

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the partitioning of the sparsity pattern columns. This is in
     * general not equal to the partitioner Epetra_Map for the domain because
     * of overlap in the matrix.
     */
    const Epetra_Map &col_partitioner () const;

    /**
     * Return a const reference to the communicator used for this object.
     */
    const Epetra_Comm &trilinos_communicator () const;
//@}
    /**
     * @name Iterators
     */
//@{

    /**
     * STL-like iterator with the first entry.
     */
    const_iterator begin () const;

    /**
     * Final iterator.
     */
    const_iterator end () const;

    /**
     * STL-like iterator with the first entry of row @p r.
     *
     * Note that if the given row is empty, i.e. does not contain any nonzero
     * entries, then the iterator returned by this function equals
     * <tt>end(r)</tt>. Note also that the iterator may not be dereferencable
     * in that case.
     */
    const_iterator begin (const size_type r) const;

    /**
     * Final iterator of row <tt>r</tt>. It points to the first element past
     * the end of line @p r, or past the end of the entire sparsity pattern.
     *
     * Note that the end iterator is not necessarily dereferencable. This is
     * in particular the case if it is the end iterator for the last row of a
     * matrix.
     */
    const_iterator end (const size_type r) const;

//@}
    /**
     * @name Input/Output
     */
//@{

    /**
     * Abstract Trilinos object that helps view in ASCII other Trilinos
     * objects. Currently this function is not implemented.  TODO: Not
     * implemented.
     */
    void write_ascii ();

    /**
     * Print (the locally owned part of) the sparsity pattern to the given
     * stream, using the format <tt>(line,col)</tt>. The optional flag outputs
     * the sparsity pattern in Trilinos style, where even the according
     * processor number is printed to the stream, as well as a summary before
     * actually writing the entries.
     */
    void print (std::ostream &out,
                const bool    write_extended_trilinos_info = false) const;

    /**
     * Print the sparsity of the matrix in a format that <tt>gnuplot</tt>
     * understands and which can be used to plot the sparsity pattern in a
     * graphical way. The format consists of pairs <tt>i j</tt> of nonzero
     * elements, each representing one entry of this matrix, one per line of
     * the output file. Indices are counted from zero on, as usual. Since
     * sparsity patterns are printed in the same way as matrices are
     * displayed, we print the negative of the column index, which means that
     * the <tt>(0,0)</tt> element is in the top left rather than in the bottom
     * left corner.
     *
     * Print the sparsity pattern in gnuplot by setting the data style to dots
     * or points and use the <tt>plot</tt> command.
     */
    void print_gnuplot (std::ostream &out) const;

//@}
    /** @addtogroup Exceptions
     * @{ */
    /**
     * Exception
     */
    DeclException1 (ExcTrilinosError,
                    int,
                    << "An error with error number " << arg1
                    << " occurred while calling a Trilinos function");

    /**
     * Exception
     */
    DeclException2 (ExcInvalidIndex,
                    size_type, size_type,
                    << "The entry with index <" << arg1 << ',' << arg2
                    << "> does not exist.");

    /**
     * Exception
     */
    DeclException0 (ExcSourceEqualsDestination);

    /**
     * Exception
     */
    DeclException4 (ExcAccessToNonLocalElement,
                    size_type, size_type, size_type, size_type,
                    << "You tried to access element (" << arg1
                    << "/" << arg2 << ")"
                    << " of a distributed matrix, but only rows "
                    << arg3 << " through " << arg4
                    << " are stored locally and can be accessed.");

    /**
     * Exception
     */
    DeclException2 (ExcAccessToNonPresentElement,
                    size_type, size_type,
                    << "You tried to access element (" << arg1
                    << "/" << arg2 << ")"
                    << " of a sparse matrix, but it appears to not"
                    << " exist in the Trilinos sparsity pattern.");
    //@}
  private:

    /**
     * Pointer to the user-supplied Epetra Trilinos mapping of the matrix
     * columns that assigns parts of the matrix to the individual processes.
     */
    std_cxx11::shared_ptr<Epetra_Map> column_space_map;

    /**
     * A sparsity pattern object in Trilinos to be used for finite element
     * based problems which allows for adding non-local elements to the
     * pattern.
     */
    std_cxx11::shared_ptr<Epetra_FECrsGraph> graph;

    /**
     * A sparsity pattern object for the non-local part of the sparsity
     * pattern that is going to be sent to the owning processor. Only used when the particular constructor or reinit method with writable_rows argument is set
     */
    std_cxx11::shared_ptr<Epetra_CrsGraph> nonlocal_graph;

    friend class SparseMatrix;
    friend class SparsityPatternIterators::Accessor;
    friend class SparsityPatternIterators::Iterator;
  };



// -------------------------- inline and template functions ----------------------


#ifndef DOXYGEN

  namespace SparsityPatternIterators
  {

    inline
    Accessor::Accessor (const SparsityPattern *sp,
                        const size_type        row,
                        const size_type        index)
      :
      sparsity_pattern(const_cast<SparsityPattern *>(sp)),
      a_row(row),
      a_index(index)
    {
      visit_present_row ();
    }


    inline
    Accessor::Accessor (const Accessor &a)
      :
      sparsity_pattern(a.sparsity_pattern),
      a_row(a.a_row),
      a_index(a.a_index),
      colnum_cache (a.colnum_cache)
    {}


    inline
    Accessor::size_type
    Accessor::row() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return a_row;
    }



    inline
    Accessor::size_type
    Accessor::column() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return (*colnum_cache)[a_index];
    }



    inline
    Accessor::size_type
    Accessor::index() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return a_index;
    }



    inline
    Iterator::Iterator(const SparsityPattern *sp,
                       const size_type        row,
                       const size_type        index)
      :
      accessor(sp, row, index)
    {}


    inline
    Iterator::Iterator(const Iterator &i)
      :
      accessor(i.accessor)
    {}



    inline
    Iterator &
    Iterator::operator++ ()
    {
      Assert (accessor.a_row < accessor.sparsity_pattern->n_rows(),
              ExcIteratorPastEnd());

      ++accessor.a_index;

      // If at end of line: do one
      // step, then cycle until we
      // find a row with a nonzero
      // number of entries.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;

          while ((accessor.a_row < accessor.sparsity_pattern->n_rows())
                 &&
                 (accessor.sparsity_pattern->row_length(accessor.a_row) == 0))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }



    inline
    Iterator
    Iterator::operator++ (int)
    {
      const Iterator old_state = *this;
      ++(*this);
      return old_state;
    }



    inline
    const Accessor &
    Iterator::operator* () const
    {
      return accessor;
    }



    inline
    const Accessor *
    Iterator::operator-> () const
    {
      return &accessor;
    }



    inline
    bool
    Iterator::operator == (const Iterator &other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    inline
    bool
    Iterator::operator != (const Iterator &other) const
    {
      return ! (*this == other);
    }



    inline
    bool
    Iterator::operator < (const Iterator &other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }

  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::begin() const
  {
    return const_iterator(this, 0, 0);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::end() const
  {
    return const_iterator(this, n_rows(), 0);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::begin(const size_type r) const
  {
    Assert (r < n_rows(), ExcIndexRangeType<size_type>(r, 0, n_rows()));
    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end (r);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::end(const size_type r) const
  {
    Assert (r < n_rows(), ExcIndexRangeType<size_type>(r, 0, n_rows()));

    // place the iterator on the first entry
    // past this line, or at the end of the
    // matrix
    for (size_type i=r+1; i<n_rows(); ++i)
      if (row_length(i) > 0)
        return const_iterator(this, i, 0);

    // if there is no such line, then take the
    // end iterator of the matrix
    return end();
  }



  inline
  bool
  SparsityPattern::in_local_range (const size_type index) const
  {
    TrilinosWrappers::types::int_type begin, end;
#ifndef DEAL_II_WITH_64BIT_INDICES
    begin = graph->RowMap().MinMyGID();
    end = graph->RowMap().MaxMyGID()+1;
#else
    begin = graph->RowMap().MinMyGID64();
    end = graph->RowMap().MaxMyGID64()+1;
#endif

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }



  inline
  bool
  SparsityPattern::is_compressed () const
  {
    return graph->Filled();
  }



  inline
  bool
  SparsityPattern::empty () const
  {
    return ((n_rows() == 0) && (n_cols() == 0));
  }



  inline
  void
  SparsityPattern::add (const size_type i,
                        const size_type j)
  {
    add_entries (i, &j, &j+1);
  }



  template <typename ForwardIterator>
  inline
  void
  SparsityPattern::add_entries (const size_type row,
                                ForwardIterator begin,
                                ForwardIterator end,
                                const bool      /*indices_are_sorted*/)
  {
    if (begin == end)
      return;

    // verify that the size of the data type Trilinos expects matches that the
    // iterator points to. we allow for some slippage between signed and
    // unsigned and only compare that they are both eiter 32 or 64 bit. to
    // write this test properly, not that we cannot compare the size of
    // '*begin' because 'begin' may be an iterator and '*begin' may be an
    // accessor class. consequently, we need to somehow get an actual value
    // from it which we can by evaluating an expression such as when
    // multiplying the value produced by 2
    Assert (sizeof(TrilinosWrappers::types::int_type) ==
            sizeof((*begin)*2),
            ExcNotImplemented());

    TrilinosWrappers::types::int_type *col_index_ptr =
      (TrilinosWrappers::types::int_type *)(&*begin);
    const int n_cols = static_cast<int>(end - begin);

    int ierr;
    if ( graph->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row)) != -1)
      ierr = graph->InsertGlobalIndices (row, n_cols, col_index_ptr);
    else if (nonlocal_graph.get() != 0)
      {
        // this is the case when we have explicitly set the off-processor rows
        // and want to create a separate matrix object for them (to retain
        // thread-safety)
        Assert (nonlocal_graph->RowMap().LID(static_cast<TrilinosWrappers::types::int_type>(row)) != -1,
                ExcMessage("Attempted to write into off-processor matrix row "
                           "that has not be specified as being writable upon "
                           "initialization"));
        ierr = nonlocal_graph->InsertGlobalIndices (row, n_cols, col_index_ptr);
      }
    else
      ierr = graph->InsertGlobalIndices
             (1, (TrilinosWrappers::types::int_type *)&row, n_cols, col_index_ptr);

    AssertThrow (ierr >= 0, ExcTrilinosError(ierr));
  }



  inline
  const Epetra_FECrsGraph &
  SparsityPattern::trilinos_sparsity_pattern () const
  {
    return *graph;
  }



  inline
  const Epetra_Map &
  SparsityPattern::domain_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->DomainMap());
  }



  inline
  const Epetra_Map &
  SparsityPattern::range_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->RangeMap());
  }



  inline
  const Epetra_Map &
  SparsityPattern::row_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->RowMap());
  }



  inline
  const Epetra_Map &
  SparsityPattern::col_partitioner () const
  {
    return static_cast<const Epetra_Map &>(graph->ColMap());
  }



  inline
  const Epetra_Comm &
  SparsityPattern::trilinos_communicator () const
  {
    return graph->RangeMap().Comm();
  }

#endif // DOXYGEN
}


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS


/*--------------------   trilinos_sparsity_pattern.h     --------------------*/

#endif
/*--------------------   trilinos_sparsity_pattern.h     --------------------*/
