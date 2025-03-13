// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_sparsity_pattern_h
#define dealii_trilinos_tpetra_sparsity_pattern_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/enable_observer_pointer.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi_stub.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/sparsity_pattern_base.h>

#  include <Tpetra_CrsGraph.hpp>

#  include <cmath>
#  include <memory>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

// forward declarations
#  ifndef DOXYGEN
class DynamicSparsityPattern;

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename MemorySpace>
    class SparsityPattern;

    template <typename Number, typename MemorySpace>
    class SparseMatrix;

    namespace SparsityPatternIterators
    {
      template <typename MemorySpace>
      class Iterator;
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    namespace SparsityPatternIterators
    {
      /**
       * Accessor class for iterators into sparsity patterns. This class is also
       * the base class for both const and non-const accessor classes into
       * sparse matrices.
       *
       * Note that this class only allows read access to elements, providing
       * their row and column number. It does not allow modifying the sparsity
       * pattern itself.
       *
       * @ingroup TrilinosWrappers
       */
      template <typename MemorySpace = dealii::MemorySpace::Host>
      class Accessor
      {
      public:
        /**
         * Declare type for container size.
         */
        using size_type = dealii::types::global_dof_index;

        /**
         * Constructor.
         */
        Accessor(const SparsityPattern<MemorySpace> *sparsity_pattern,
                 const size_type                     row,
                 const size_type                     index);

        /**
         * Row number of the element represented by this object.
         */
        size_type
        row() const;

        /**
         * Index in row of the element represented by this object.
         */
        size_type
        index() const;

        /**
         * Column number of the element represented by this object.
         */
        size_type
        column() const;

        /**
         * Exception
         */
        DeclException0(ExcBeyondEndOfSparsityPattern);

        /**
         * Exception
         */
        DeclException3(ExcAccessToNonlocalRow,
                       size_type,
                       size_type,
                       size_type,
                       << "You tried to access row " << arg1
                       << " of a distributed sparsity pattern, "
                       << " but only rows " << arg2 << " through " << arg3
                       << " are stored locally and can be accessed.");

      private:
        /**
         * The matrix accessed.
         */
        SparsityPattern<MemorySpace> *sparsity_pattern;

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
        std::shared_ptr<std::vector<dealii::types::signed_global_dof_index>>
          colnum_cache;

        /**
         * Discard the old row caches (they may still be used by other
         * accessors) and generate new ones for the row pointed to presently by
         * this accessor.
         */
        void
        visit_present_row();

        // Make enclosing class a friend.
        friend class Iterator<MemorySpace>;
      };

      /**
       * Iterator class for sparsity patterns of type
       * TrilinosWrappers::SparsityPattern. Access to individual elements of the
       * sparsity pattern is handled by the Accessor class in this namespace.
       */
      template <typename MemorySpace = dealii::MemorySpace::Host>
      class Iterator
      {
      public:
        /**
         * Declare type for container size.
         */
        using size_type = types::global_dof_index;

        /**
         * Constructor. Create an iterator into the matrix @p matrix for the
         * given row and the index within it.
         */
        Iterator(const SparsityPattern<MemorySpace> *sparsity_pattern,
                 const size_type                     row,
                 const size_type                     index);

        /**
         * Copy constructor.
         */
        Iterator(const Iterator<MemorySpace> &i);

        /**
         * Prefix increment.
         */
        Iterator<MemorySpace> &
        operator++();

        /**
         * Postfix increment.
         */
        Iterator
        operator++(int);

        /**
         * Dereferencing operator.
         */
        const Accessor<MemorySpace> &
        operator*() const;

        /**
         * Dereferencing operator.
         */
        const Accessor<MemorySpace> *
        operator->() const;

        /**
         * Comparison. True, if both iterators point to the same matrix
         * position.
         */
        bool
        operator==(const Iterator<MemorySpace> &) const;

        /**
         * Inverse of <tt>==</tt>.
         */
        bool
        operator!=(const Iterator<MemorySpace> &) const;

        /**
         * Comparison operator. Result is true if either the first row number is
         * smaller or if the row numbers are equal and the first index is
         * smaller.
         */
        bool
        operator<(const Iterator<MemorySpace> &) const;

        /**
         * Exception
         */
        DeclException2(ExcInvalidIndexWithinRow,
                       size_type,
                       size_type,
                       << "Attempt to access element " << arg2 << " of row "
                       << arg1 << " which doesn't have that many elements.");

      private:
        /**
         * Store an object of the accessor class.
         */
        Accessor<MemorySpace> accessor;

        friend class TpetraWrappers::SparsityPattern<MemorySpace>;
      };

    } // namespace SparsityPatternIterators


    /**
     * This class implements a wrapper class to use the Trilinos distributed
     * sparsity pattern class Tpetra::CrsGraph. This class is designed to be
     * used for construction of %parallel Trilinos matrices. The functionality
     * of this class is modeled after the existing sparsity pattern classes,
     * with the difference that this class can work fully in %parallel according
     * to a partitioning of the sparsity pattern rows.
     *
     * This class has many similarities to the  DynamicSparsityPattern, since it
     * can dynamically add elements to the pattern without any memory being
     * previously reserved for it. However, it also has a method
     * SparsityPattern<MemorySpace>::compress(), that finalizes the pattern and
     * enables its use with Trilinos sparse matrices.
     *
     * @ingroup TrilinosWrappers
     * @ingroup Sparsity
     */
    template <typename MemorySpace = dealii::MemorySpace::Host>
    class SparsityPattern : public SparsityPatternBase
    {
    public:
      /**
       * Declare type for container size.
       */
      using size_type = dealii::types::global_dof_index;
      /**
       * Declare an alias for the iterator class.
       */
      using const_iterator = SparsityPatternIterators::Iterator<MemorySpace>;

      /**
       * @name Basic constructors and initialization
       */
      /** @{ */
      /**
       * Default constructor. Generates an empty (zero-size) sparsity pattern.
       */
      SparsityPattern();

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
       * for the deal.II DynamicSparsityPattern classes. However, a good
       * estimate will reduce the setup time of the sparsity pattern.
       */
      SparsityPattern(const size_type m,
                      const size_type n,
                      const size_type n_entries_per_row = 0);

      /**
       * Generate a sparsity pattern that is completely stored locally, having
       * $m$ rows and $n$ columns. The resulting matrix will be completely
       * stored locally, too.
       *
       * The vector <tt>n_entries_per_row</tt> specifies the number of entries
       * in each row (an information usually not available, though).
       */
      SparsityPattern(const size_type               m,
                      const size_type               n,
                      const std::vector<size_type> &n_entries_per_row);

      /**
       * Move constructor. Create a new sparse matrix by stealing the internal
       * data.
       */
      SparsityPattern(SparsityPattern<MemorySpace> &&other) noexcept;

      /**
       * Copy constructor. Sets the calling sparsity pattern to be the same as
       * the input sparsity pattern.
       */
      SparsityPattern(
        const SparsityPattern<MemorySpace> &input_sparsity_pattern);

      /**
       * Destructor. Made virtual so that one can use pointers to this class.
       */
      virtual ~SparsityPattern() override = default;

      /**
       * Initialize a sparsity pattern that is completely stored locally, having
       * $m$ rows and $n$ columns. The resulting matrix will be completely
       * stored locally.
       *
       * The number of columns entries per row is specified as the maximum
       * number of entries argument.  This does not need to be an accurate
       * number since the entries are allocated dynamically in a similar manner
       * as for the deal.II DynamicSparsityPattern classes, but a good estimate
       * will reduce the setup time of the sparsity pattern.
       */
      void
      reinit(const size_type m,
             const size_type n,
             const size_type n_entries_per_row = 0);

      /**
       * Initialize a sparsity pattern that is completely stored locally, having
       * $m$ rows and $n$ columns. The resulting matrix will be completely
       * stored locally.
       *
       * The vector <tt>n_entries_per_row</tt> specifies the number of entries
       * in each row.
       */
      void
      reinit(const size_type               m,
             const size_type               n,
             const std::vector<size_type> &n_entries_per_row);

      /**
       * Copy function. Sets the calling sparsity pattern to be the same as the
       * input sparsity pattern.
       */
      void
      copy_from(const SparsityPattern<MemorySpace> &input_sparsity_pattern);

      /**
       * Copy function from one of the deal.II sparsity patterns. If used in
       * parallel, this function uses an ad-hoc partitioning of the rows and
       * columns.
       */
      template <typename SparsityPatternType>
      void
      copy_from(const SparsityPatternType &nontrilinos_sparsity_pattern);

      /**
       * Copy operator. This operation is only allowed for empty objects, to
       * avoid potentially very costly operations automatically synthesized by
       * the compiler. Use copy_from() instead if you know that you really want
       * to copy a sparsity pattern with non-trivial content.
       */
      SparsityPattern<MemorySpace> &
      operator=(const SparsityPattern<MemorySpace> &input_sparsity_pattern);

      /**
       * Release all memory and return to a state just like after having called
       * the default constructor.
       *
       * This is a collective operation that needs to be called on all
       * processors in order to avoid a dead lock.
       */
      void
      clear();

      /**
       * In analogy to our own SparsityPattern class, this function compresses
       * the sparsity pattern and allows the resulting pattern to be used for
       * actually generating a (Trilinos-based) matrix. This function also
       * exchanges non-local data that might have accumulated during the
       * addition of new elements. This function must therefore be called once
       * the structure is fixed. This is a collective operation, i.e., it needs
       * to be run on all processors when used in parallel.
       */
      void
      compress();
      /** @} */

      /**
       * @name Constructors and initialization using an IndexSet description
       */
      /** @{ */

      /**
       * Constructor for a square sparsity pattern using an IndexSet and an MPI
       * communicator for the description of the %parallel partitioning.
       * Moreover, the number of nonzero entries in the rows of the sparsity
       * pattern can be specified. Note that this number does not need to be
       * exact, and it is even allowed that the actual sparsity structure has
       * more nonzero entries than specified in the constructor. However it is
       * still advantageous to provide good estimates here since a good value
       * will avoid repeated allocation of memory, which considerably increases
       * the performance when creating the sparsity pattern.
       */
      SparsityPattern(const IndexSet &parallel_partitioning,
                      const MPI_Comm  communicator      = MPI_COMM_WORLD,
                      const size_type n_entries_per_row = 0);

      /**
       * Same as before, but now use the exact number of nonzero entries in
       * each m row. Since we know the number of elements in the sparsity
       * pattern exactly in this case, we can already allocate the right amount
       * of memory, which makes the creation process by the respective
       * SparsityPattern<MemorySpace>::reinit call considerably faster. However,
       * this is a rather unusual situation, since knowing the number of entries
       * in each row is usually connected to knowing the indices of nonzero
       * entries, which the sparsity pattern is designed to describe.
       */
      SparsityPattern(const IndexSet               &parallel_partitioning,
                      const MPI_Comm                communicator,
                      const std::vector<size_type> &n_entries_per_row);

      /**
       * This constructor is similar to the one above, but it now takes two
       * different index sets to describe the %parallel partitioning of rows and
       * columns. This interface is meant to be used for generating rectangular
       * sparsity pattern. Note that there is no real parallelism along the
       * columns &ndash; the processor that owns a certain row always owns all
       * the column elements, no matter how far they might be spread out. The
       * second Tpetra::Map is only used to specify the number of columns and
       * for internal arrangements when doing matrix-vector products with
       * vectors based on that column map.
       *
       * The number of columns entries per row is specified as the maximum
       * number of entries argument.
       */
      SparsityPattern(const IndexSet &row_parallel_partitioning,
                      const IndexSet &col_parallel_partitioning,
                      const MPI_Comm  communicator      = MPI_COMM_WORLD,
                      const size_type n_entries_per_row = 0);

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
      SparsityPattern(const IndexSet               &row_parallel_partitioning,
                      const IndexSet               &col_parallel_partitioning,
                      const MPI_Comm                communicator,
                      const std::vector<size_type> &n_entries_per_row);

      /**
       * This constructor constructs general sparsity patterns, possible non-
       * square ones. Constructing a sparsity pattern this way allows the user
       * to explicitly specify the rows into which we are going to add elements.
       * This set is required to be a superset of the first index set @p
       * row_parallel_partitioning that includes also rows that are owned by
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
       * this flexibility come at a cost, the most prominent being that adding
       * elements into the same matrix from multiple threads in shared memory is
       * not safe whenever MPI is used. For these settings, the current method
       * is the one to choose: It will store the off-processor data as an
       * additional sparsity pattern (that is then passed to the Trilinos matrix
       * via the reinit method) which can be organized in such a way that
       * thread-safety can be ensured (as long as the user makes sure to never
       * write into the same matrix row simultaneously, of course).
       */
      SparsityPattern(const IndexSet &row_parallel_partitioning,
                      const IndexSet &col_parallel_partitioning,
                      const IndexSet &writable_rows,
                      const MPI_Comm  communicator      = MPI_COMM_WORLD,
                      const size_type n_entries_per_row = 0);

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
      reinit(const IndexSet &parallel_partitioning,
             const MPI_Comm  communicator      = MPI_COMM_WORLD,
             const size_type n_entries_per_row = 0);

      /**
       * Same as before, but now use the exact number of nonzero entries in
       * each row. Since we know the number of elements in the sparsity pattern
       * exactly in this case, we can already allocate the right amount of
       * memory, which makes process of adding entries to the sparsity pattern
       * considerably faster. However, this is a rather unusual situation, since
       * knowing the number of entries in each row is usually connected to
       * knowing the indices of nonzero entries, which the sparsity pattern is
       * designed to describe.
       */
      void
      reinit(const IndexSet               &parallel_partitioning,
             const MPI_Comm                communicator,
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
       * arrangements when doing matrix-vector products with vectors based on an
       * Tpetra::Map based on that IndexSet.
       *
       * The number of columns entries per row is specified by the argument
       * <tt>n_entries_per_row</tt>.
       */
      void
      reinit(const IndexSet &row_parallel_partitioning,
             const IndexSet &col_parallel_partitioning,
             const MPI_Comm  communicator      = MPI_COMM_WORLD,
             const size_type n_entries_per_row = 0);

      /**
       * This reinit function is used to specify general matrices, possibly non-
       * square ones. In addition to the arguments of the other reinit method
       * above, it allows the user to explicitly specify the rows into which we
       * are going to add elements. This set is a superset of the first index
       * set @p row_parallel_partitioning that includes also rows that are owned
       * by another processor (ghost rows).
       *
       * This method is beneficial when the rows to which a processor is going
       * to write can be determined before actually inserting elements into the
       * matrix. For the typical parallel::distributed::Triangulation class used
       * in deal.II, we know that a processor only will add row elements for
       * what we call the locally relevant dofs (see
       * DoFTools::extract_locally_relevant_dofs). Trilinos matrices allow to
       * add elements to arbitrary rows (as done by all the other reinit
       * functions) and this is what all the other reinit methods do, too.
       * However, this flexibility come at a cost, the most prominent being that
       * adding elements into the same matrix from multiple threads in shared
       * memory is not safe whenever MPI is used. For these settings, the
       * current method is the one to choose: It will store the off-processor
       * data as an additional sparsity pattern (that is then passed to the
       * Trilinos matrix via the reinit method) which can be organized in such a
       * way that thread-safety can be ensured (as long as the user makes sure
       * to never write into the same matrix row simultaneously, of course).
       */
      void
      reinit(const IndexSet &row_parallel_partitioning,
             const IndexSet &col_parallel_partitioning,
             const IndexSet &writeable_rows,
             const MPI_Comm  communicator      = MPI_COMM_WORLD,
             const size_type n_entries_per_row = 0);

      /**
       * Same as before, but now using a vector <tt>n_entries_per_row</tt> for
       * specifying the number of entries in each row of the sparsity pattern.
       */
      void
      reinit(const IndexSet               &row_parallel_partitioning,
             const IndexSet               &col_parallel_partitioning,
             const MPI_Comm                communicator,
             const std::vector<size_type> &n_entries_per_row);

      /**
       * Reinit function. Takes one of the deal.II sparsity patterns and the
       * %parallel partitioning of the rows and columns specified by two index
       * sets and a %parallel communicator for initializing the current Trilinos
       * sparsity pattern. The optional argument @p exchange_data can be used
       * for reinitialization with a sparsity pattern that is not fully
       * constructed. This feature is only implemented for input sparsity
       * patterns of type DynamicSparsityPattern.
       */
      template <typename SparsityPatternType>
      void
      reinit(const IndexSet            &row_parallel_partitioning,
             const IndexSet            &col_parallel_partitioning,
             const SparsityPatternType &nontrilinos_sparsity_pattern,
             const MPI_Comm             communicator  = MPI_COMM_WORLD,
             const bool                 exchange_data = false);

      /**
       * Reinit function. Takes one of the deal.II sparsity patterns and a
       * %parallel partitioning of the rows and columns for initializing the
       * current Trilinos sparsity pattern. The optional argument @p
       * exchange_data can be used for reinitialization with a sparsity pattern
       * that is not fully constructed. This feature is only implemented for
       * input sparsity patterns of type DynamicSparsityPattern.
       */
      template <typename SparsityPatternType>
      void
      reinit(const IndexSet            &parallel_partitioning,
             const SparsityPatternType &nontrilinos_sparsity_pattern,
             const MPI_Comm             communicator  = MPI_COMM_WORLD,
             const bool                 exchange_data = false);
      /** @} */
      /**
       * @name Information on the sparsity pattern
       */
      /** @{ */

      /**
       * Return the state of the sparsity pattern, i.e., whether compress()
       * needs to be called after an operation requiring data exchange.
       */
      bool
      is_compressed() const;

      /**
       * Return the maximum number of entries per row on the current processor.
       */
      unsigned int
      max_entries_per_row() const;

      /**
       * Return the local dimension of the sparsity pattern, i.e. the number of
       * rows stored on the present MPI process. In the sequential case, this
       * number is the same as n_rows(), but for parallel matrices it may be
       * smaller.
       *
       * To figure out which elements exactly are stored locally, use
       * local_range().
       */
      unsigned int
      local_size() const;

      /**
       * Return a pair of indices indicating which rows of this sparsity pattern
       * are stored locally. The first number is the index of the first row
       * stored, the second the index of the one past the last one that is
       * stored locally. If this is a sequential matrix, then the result will be
       * the pair (0,n_rows()), otherwise it will be a pair (i,i+n), where
       * <tt>n=local_size()</tt>.
       */
      std::pair<size_type, size_type>
      local_range() const;

      /**
       * Return whether @p index is in the local range or not, see also
       * local_range().
       */
      bool
      in_local_range(const size_type index) const;

      /**
       * Return the number of nonzero elements of this sparsity pattern.
       */
      std::uint64_t
      n_nonzero_elements() const;

      /**
       * Return the number of entries in the given row.
       *
       * In a parallel context, the row in question may of course not be
       * stored on the current processor, and in that case it is not
       * possible to query the number of entries in it. In that case,
       * the returned value is `static_cast<size_type>(-1)`.
       */
      size_type
      row_length(const size_type row) const;

      /**
       * Compute the bandwidth of the matrix represented by this structure. The
       * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
       * represents a nonzero entry of the matrix. Consequently, the maximum
       * bandwidth a $n\times m$ matrix can have is $\max\{n-1,m-1\}$.
       */
      size_type
      bandwidth() const;

      /**
       * Return whether the object is empty. It is empty if no memory is
       * allocated, which is the same as when both dimensions are zero.
       */
      bool
      empty() const;

      /**
       * Return whether the index (<i>i,j</i>) exists in the sparsity pattern
       * (i.e., it may be nonzero) or not.
       */
      bool
      exists(const size_type i, const size_type j) const;

      /**
       * Return whether a given @p row is stored in the current object
       * on this process.
       */
      bool
      row_is_stored_locally(const size_type i) const;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object. Currently not implemented for this class.
       */
      std::size_t
      memory_consumption() const;

      /** @} */
      /**
       * @name Adding entries
       */
      /** @{ */
      /**
       * Add the element (<i>i,j</i>) to the sparsity pattern.
       */
      void
      add(const size_type i, const size_type j);


      /**
       * Add several elements in one row to the sparsity pattern.
       */
      template <typename ForwardIterator>
      void
      add_entries(const size_type row,
                  ForwardIterator begin,
                  ForwardIterator end,
                  const bool      indices_are_sorted = false);

      virtual void
      add_row_entries(
        const dealii::types::global_dof_index                  &row,
        const ArrayView<const dealii::types::global_dof_index> &columns,
        const bool indices_are_sorted = false) override;

      using SparsityPatternBase::add_entries;

      /** @} */
      /**
       * @name Access of underlying Trilinos data
       */
      /** @{ */

      /**
       * Return a Teuchos::RCP to the underlying Trilinos Tpetra::CrsGraph
       * data that stores the sparsity pattern.
       */
      Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>>
      trilinos_sparsity_pattern() const;

      /**
       * Return a const Teuchos::RCP to the underlying Trilinos Tpetra::Map that
       * sets the parallel partitioning of the domain space of this sparsity
       * pattern, i.e., the partitioning of the vectors matrices based on this
       * sparsity pattern are multiplied with.
       */
      Teuchos::RCP<const TpetraTypes::MapType<MemorySpace>>
      domain_partitioner() const;

      /**
       * Return a const Teuchos::RCP to the underlying Trilinos Tpetra::Map that
       * sets the partitioning of the range space of this sparsity pattern,
       * i.e., the partitioning of the vectors that are result from matrix-
       * vector products.
       */
      Teuchos::RCP<const TpetraTypes::MapType<MemorySpace>>
      range_partitioner() const;

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;

      /**
       * Return the underlying Teuchos::MPI communicator.
       */
      Teuchos::RCP<const Teuchos::Comm<int>>
      get_teuchos_mpi_communicator() const;

      /** @} */

      /**
       * @name Partitioners
       */
      /** @{ */

      /**
       * Return the partitioning of the domain space of this pattern, i.e., the
       * partitioning of the vectors a matrix based on this sparsity pattern has
       * to be multiplied with.
       */
      IndexSet
      locally_owned_domain_indices() const;

      /**
       * Return the partitioning of the range space of this pattern, i.e., the
       * partitioning of the vectors that are the result from matrix-vector
       * products from a matrix based on this pattern.
       */
      IndexSet
      locally_owned_range_indices() const;

      /** @} */

      /**
       * @name Iterators
       */
      /** @{ */

      /**
       * Iterator starting at the first entry.
       */
      const_iterator
      begin() const;

      /**
       * Final iterator.
       */
      const_iterator
      end() const;

      /**
       * Iterator starting at the first entry of row @p r.
       *
       * Note that if the given row is empty, i.e. does not contain any nonzero
       * entries, then the iterator returned by this function equals
       * <tt>end(r)</tt>. Note also that the iterator may not be dereferenceable
       * in that case.
       */
      const_iterator
      begin(const size_type r) const;

      /**
       * Final iterator of row <tt>r</tt>. It points to the first element past
       * the end of line @p r, or past the end of the entire sparsity pattern.
       *
       * Note that the end iterator is not necessarily dereferenceable. This is
       * in particular the case if it is the end iterator for the last row of a
       * matrix.
       */
      const_iterator
      end(const size_type r) const;

      /** @} */
      /**
       * @name Input/Output
       */
      /** @{ */

      /**
       * Print (the locally owned part of) the sparsity pattern to the given
       * stream, using the format <tt>(line,col)</tt>. The optional flag outputs
       * the sparsity pattern in Trilinos style, where even the according
       * processor number is printed to the stream, as well as a summary before
       * actually writing the entries.
       */
      void
      print(std::ostream &out,
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
      void
      print_gnuplot(std::ostream &out) const;

      /** @} */
      /**
       * @addtogroup Exceptions
       * @{
       */
      /**
       * Exception
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

      /**
       * Exception
       */
      DeclException2(ExcInvalidIndex,
                     size_type,
                     size_type,
                     << "The entry with index <" << arg1 << ',' << arg2
                     << "> does not exist.");

      /**
       * Exception
       */
      DeclException4(ExcAccessToNonLocalElement,
                     size_type,
                     size_type,
                     size_type,
                     size_type,
                     << "You tried to access element (" << arg1 << '/' << arg2
                     << ')'
                     << " of a distributed matrix, but only rows in range ["
                     << arg3 << ',' << arg4
                     << "] are stored locally and can be accessed.");

      /**
       * Exception
       */
      DeclException2(ExcAccessToNonPresentElement,
                     size_type,
                     size_type,
                     << "You tried to access element (" << arg1 << '/' << arg2
                     << ')' << " of a sparse matrix, but it appears to not"
                     << " exist in the Trilinos sparsity pattern.");
      /** @} */
    private:
      /**
       * Teuchos::RCP to the user-supplied Tpetra Trilinos mapping of the matrix
       * columns that assigns parts of the matrix to the individual processes.
       */
      Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> column_space_map;

      /**
       * A sparsity pattern object in Trilinos to be used for finite element
       * based problems which allows for adding non-local elements to the
       * pattern.
       */
      Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> graph;

      /**
       * A sparsity pattern object for the non-local part of the sparsity
       * pattern that is going to be sent to the owning processor. Only used
       * when the particular constructor or reinit method with writable_rows
       * argument is set
       */
      Teuchos::RCP<TpetraTypes::GraphType<MemorySpace>> nonlocal_graph;

      // TODO: currently only for double
      friend class SparseMatrix<double, MemorySpace>;
      friend class SparsityPatternIterators::Accessor<MemorySpace>;
      friend class SparsityPatternIterators::Iterator<MemorySpace>;
    };



    // ---------------- inline and template functions -----------------


#  ifndef DOXYGEN

    namespace SparsityPatternIterators
    {
      template <typename MemorySpace>
      inline Accessor<MemorySpace>::Accessor(
        const SparsityPattern<MemorySpace> *sp,
        const size_type                     row,
        const size_type                     index)
        : sparsity_pattern(const_cast<SparsityPattern<MemorySpace> *>(sp))
        , a_row(row)
        , a_index(index)
      {
        visit_present_row();
      }



      template <typename MemorySpace>
      inline typename Accessor<MemorySpace>::size_type
      Accessor<MemorySpace>::row() const
      {
        Assert(a_row < sparsity_pattern->n_rows(),
               ExcBeyondEndOfSparsityPattern());
        return a_row;
      }



      template <typename MemorySpace>
      inline typename Accessor<MemorySpace>::size_type
      Accessor<MemorySpace>::column() const
      {
        Assert(a_row < sparsity_pattern->n_rows(),
               ExcBeyondEndOfSparsityPattern());
        return (*colnum_cache)[a_index];
      }



      template <typename MemorySpace>
      inline typename Accessor<MemorySpace>::size_type
      Accessor<MemorySpace>::index() const
      {
        Assert(a_row < sparsity_pattern->n_rows(),
               ExcBeyondEndOfSparsityPattern());
        return a_index;
      }



      template <typename MemorySpace>
      inline Iterator<MemorySpace>::Iterator(
        const SparsityPattern<MemorySpace> *sp,
        const size_type                     row,
        const size_type                     index)
        : accessor(sp, row, index)
      {}



      template <typename MemorySpace>
      inline Iterator<MemorySpace>::Iterator(const Iterator<MemorySpace> &) =
        default;



      template <typename MemorySpace>
      inline Iterator<MemorySpace> &
      Iterator<MemorySpace>::operator++()
      {
        Assert(accessor.a_row < accessor.sparsity_pattern->n_rows(),
               ExcIteratorPastEnd());

        ++accessor.a_index;

        // If at end of line: do one step, then cycle until we find a row with a
        // nonzero number of entries that is stored locally.
        if (accessor.a_index >=
            static_cast<dealii::types::signed_global_dof_index>(
              accessor.colnum_cache->size()))
          {
            accessor.a_index = 0;
            ++accessor.a_row;

            while (accessor.a_row <
                   static_cast<dealii::types::signed_global_dof_index>(
                     accessor.sparsity_pattern->n_rows()))
              {
                const auto row_length =
                  accessor.sparsity_pattern->row_length(accessor.a_row);
                if (row_length == 0 ||
                    !accessor.sparsity_pattern->row_is_stored_locally(
                      accessor.a_row))
                  ++accessor.a_row;
                else
                  break;
              }

            accessor.visit_present_row();
          }
        return *this;
      }



      template <typename MemorySpace>
      inline Iterator<MemorySpace>
      Iterator<MemorySpace>::operator++(int)
      {
        const Iterator<MemorySpace> old_state = *this;
        ++(*this);
        return old_state;
      }



      template <typename MemorySpace>
      inline const Accessor<MemorySpace> &
      Iterator<MemorySpace>::operator*() const
      {
        return accessor;
      }



      template <typename MemorySpace>
      inline const Accessor<MemorySpace> *
      Iterator<MemorySpace>::operator->() const
      {
        return &accessor;
      }



      template <typename MemorySpace>
      inline bool
      Iterator<MemorySpace>::operator==(
        const Iterator<MemorySpace> &other) const
      {
        return (accessor.a_row == other.accessor.a_row &&
                accessor.a_index == other.accessor.a_index);
      }



      template <typename MemorySpace>
      inline bool
      Iterator<MemorySpace>::operator!=(
        const Iterator<MemorySpace> &other) const
      {
        return !(*this == other);
      }



      template <typename MemorySpace>
      inline bool
      Iterator<MemorySpace>::operator<(const Iterator<MemorySpace> &other) const
      {
        return (accessor.row() < other.accessor.row() ||
                (accessor.row() == other.accessor.row() &&
                 accessor.index() < other.accessor.index()));
      }

    } // namespace SparsityPatternIterators



    template <typename MemorySpace>
    inline typename SparsityPattern<MemorySpace>::const_iterator
    SparsityPattern<MemorySpace>::begin() const
    {
      const size_type first_valid_row = this->local_range().first;
      return const_iterator(this, first_valid_row, 0);
    }



    template <typename MemorySpace>
    inline typename SparsityPattern<MemorySpace>::const_iterator
    SparsityPattern<MemorySpace>::end() const
    {
      return const_iterator(this, n_rows(), 0);
    }



    template <typename MemorySpace>
    inline typename SparsityPattern<MemorySpace>::const_iterator
    SparsityPattern<MemorySpace>::begin(const size_type r) const
    {
      AssertIndexRange(r, n_rows());
      if (row_length(r) > 0)
        return const_iterator(this, r, 0);
      else
        return end(r);
    }



    template <typename MemorySpace>
    inline typename SparsityPattern<MemorySpace>::const_iterator
    SparsityPattern<MemorySpace>::end(const size_type r) const
    {
      AssertIndexRange(r, n_rows());

      // place the iterator on the first entry
      // past this line, or at the end of the
      // matrix
      for (size_type i = r + 1; i < n_rows(); ++i)
        if (row_length(i) > 0)
          return const_iterator(this, i, 0);

      // if there is no such line, then take the
      // end iterator of the matrix
      return end();
    }



    template <typename MemorySpace>
    inline bool
    SparsityPattern<MemorySpace>::in_local_range(const size_type index) const
    {
      const TrilinosWrappers::types::int_type begin =
        graph->getRowMap()->getMinGlobalIndex();
      const TrilinosWrappers::types::int_type end =
        graph->getRowMap()->getMaxGlobalIndex() + 1;

      return ((index >= static_cast<size_type>(begin)) &&
              (index < static_cast<size_type>(end)));
    }



    template <typename MemorySpace>
    inline bool
    SparsityPattern<MemorySpace>::is_compressed() const
    {
      return graph->isFillComplete();
    }



    template <typename MemorySpace>
    inline bool
    SparsityPattern<MemorySpace>::empty() const
    {
      return ((n_rows() == 0) && (n_cols() == 0));
    }



    template <typename MemorySpace>
    inline void
    SparsityPattern<MemorySpace>::add(const size_type i, const size_type j)
    {
      add_entries(i, &j, &j + 1);
    }



    template <typename MemorySpace>
    template <typename ForwardIterator>
    inline void
    SparsityPattern<MemorySpace>::add_entries(const size_type row,
                                              ForwardIterator begin,
                                              ForwardIterator end,
                                              const bool /*indices_are_sorted*/)
    {
      if (begin == end)
        return;

      // verify that the size of the data type Trilinos expects matches that the
      // iterator points to. we allow for some slippage between signed and
      // unsigned and only compare that they are both either 32 or 64 bit. to
      // write this test properly, not that we cannot compare the size of
      // '*begin' because 'begin' may be an iterator and '*begin' may be an
      // accessor class. consequently, we need to somehow get an actual value
      // from it which we can by evaluating an expression such as when
      // multiplying the value produced by 2
      Assert(sizeof(TrilinosWrappers::types::int_type) == sizeof((*begin) * 2),
             ExcNotImplemented());

      const TrilinosWrappers::types::int_type *col_index_ptr_begin =
        reinterpret_cast<TrilinosWrappers::types::int_type *>(
          const_cast<std::decay_t<decltype(*begin)> *>(&*begin));

      const TrilinosWrappers::types::int_type *col_index_ptr_end =
        reinterpret_cast<TrilinosWrappers::types::int_type *>(
          const_cast<std::decay_t<decltype(*end)> *>(&*end));

      // Check at least for the first index that the conversion actually works
      AssertDimension(*col_index_ptr_begin, *begin);
      AssertDimension(*col_index_ptr_end, *end);
      TrilinosWrappers::types::int_type trilinos_row_index = row;

      // TODO: The following line creates an array by copying the entries.
      //       Perhaps there is a way to only create a 'view' of these arrays
      //       and pass that to Tpetra?
      Teuchos::Array<TrilinosWrappers::types::int_type> array(
        col_index_ptr_begin, col_index_ptr_end);

      if (row_is_stored_locally(row))
        graph->insertGlobalIndices(trilinos_row_index, array());
      else if (nonlocal_graph.get() != nullptr)
        {
          // this is the case when we have explicitly set the off-processor rows
          // and want to create a separate matrix object for them (to retain
          // thread-safety)
          Assert(nonlocal_graph->getRowMap()->getLocalElement(row) !=
                   Teuchos::OrdinalTraits<
                     dealii::types::signed_global_dof_index>::invalid(),
                 ExcMessage("Attempted to write into off-processor matrix row "
                            "that has not be specified as being writable upon "
                            "initialization"));
          nonlocal_graph->insertGlobalIndices(trilinos_row_index, array);
        }
      else
        graph->insertGlobalIndices(trilinos_row_index, array);
    }



    template <typename MemorySpace>
    inline Teuchos::RCP<Tpetra::CrsGraph<int,
                                         dealii::types::signed_global_dof_index,
                                         TpetraTypes::NodeType<MemorySpace>>>
    SparsityPattern<MemorySpace>::trilinos_sparsity_pattern() const
    {
      return graph;
    }



    template <typename MemorySpace>
    inline IndexSet
    SparsityPattern<MemorySpace>::locally_owned_domain_indices() const
    {
      return IndexSet(graph->getDomainMap().getConst());
    }



    template <typename MemorySpace>
    inline IndexSet
    SparsityPattern<MemorySpace>::locally_owned_range_indices() const
    {
      return IndexSet(graph->getRangeMap().getConst());
    }

#  endif // DOXYGEN
  }      // namespace TpetraWrappers

} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif
