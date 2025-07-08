// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_read_write_vector_h
#define dealii_read_write_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/read_vector.h>
#include <deal.II/lac/vector_operation.h>

#include <cstdlib>
#include <cstring>
#include <iomanip>

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/trilinos_epetra_vector.h>
#  include <deal.II/lac/trilinos_tpetra_block_vector.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_MultiVector.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class Vector;

namespace LinearAlgebra
{
  template <typename>
  class Vector;
  namespace distributed
  {
    template <typename, typename>
    class Vector;
  } // namespace distributed
} // namespace LinearAlgebra

#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace TrilinosWrappers
#  endif

#endif

namespace LinearAlgebra
{
  /**
   * @addtogroup Vectors
   * @{
   */

  /**
   * ReadWriteVector is intended to represent vectors in ${\mathbb R}^N$ for
   * which it stores all or a subset of elements. The latter case in important
   * in parallel computations, where $N$ may be so large that no processor can
   * actually store all elements of a solution vector, but where this is also
   * not necessary: one typically only has to store the values of degrees of
   * freedom that live on cells that are locally owned plus potentially those
   * degrees of freedom that live on ghost cells.
   *
   * This class provides access to individual elements to be read or written.
   * However, it does not allow global operations such as taking the norm
   * or dot products between vectors.
   *
   * <h3>Storing elements</h3>
   * Most of the time, one will simply read from or
   * write into a vector of the current class using the global numbers of
   * these degrees of freedom. This is done using operator()() or operator[]()
   * which call global_to_local() to transform the <i>global</i> index into a
   * <i>local</i> one. In such cases, it is clear that one can only access
   * elements of the vector that the current object indeed stores.
   *
   * However, it is also possible to access elements in the order in which
   * they are stored by the current object. In other words, one is not
   * interested in accessing elements with their <i>global</i> indices, but
   * instead using an enumeration that only takes into account the elements
   * that are actually stored. This is facilitated by the local_element()
   * function. To this end, it is necessary to know <i>in which order</i> the
   * current class stores its element. The elements of all the consecutive
   * ranges are stored in ascending order of the first index of each range.
   * The function IndexSet::largest_range_starting_index() can be used to
   * get the first index of the largest range.
   */
  template <typename Number>
  class ReadWriteVector : public ReadVector<Number>
  {
  public:
    /**
     * Declare standard types used in all containers. These types parallel
     * those in the <tt>C++</tt> standard libraries <tt>vector<...></tt>
     * class.
     */
    using value_type      = Number;
    using pointer         = value_type *;
    using const_pointer   = const value_type *;
    using iterator        = value_type *;
    using const_iterator  = const value_type *;
    using reference       = value_type &;
    using const_reference = const value_type &;
    using size_type       = types::global_dof_index;
    using real_type       = typename numbers::NumberTraits<Number>::real_type;

    /**
     * @name 1: Basic Object-handling
     */
    /** @{ */
    /**
     * Empty constructor.
     */
    ReadWriteVector();

    /**
     * Copy constructor.
     */
    ReadWriteVector(const ReadWriteVector<Number> &in_vector);

    /**
     * Construct a vector given the size, the stored elements have their
     * index in [0,size).
     */
    explicit ReadWriteVector(const size_type size);

    /**
     * Construct a vector whose stored elements indices are given by the
     * IndexSet @p locally_stored_indices.
     */
    explicit ReadWriteVector(const IndexSet &locally_stored_indices);

    /**
     * Destructor.
     */
    ~ReadWriteVector() override = default;

    /**
     * Set the global size of the vector to @p size. The stored elements have
     * their index in [0,size).
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    virtual void
    reinit(const size_type size, const bool omit_zeroing_entries = false);

    /**
     * Uses the same IndexSet as the one of the input vector @p in_vector and
     * allocates memory for this vector.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    template <typename Number2>
    void
    reinit(const ReadWriteVector<Number2> &in_vector,
           const bool                      omit_zeroing_entries = false);

    /**
     * Initializes the vector. The indices are specified by @p
     * locally_stored_indices.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     * locally_stored_indices.
     */
    virtual void
    reinit(const IndexSet &locally_stored_indices,
           const bool      omit_zeroing_entries = false);


#ifdef DEAL_II_WITH_TRILINOS
    /**
     * Initialize this ReadWriteVector by supplying access to all locally
     * available entries in the given ghosted or non-ghosted vector.
     *
     * @note This function currently copies the values from the argument into
     * the ReadWriteVector, so modifications here will not modify @p trilinos_vec.
     *
     * This function is mainly written for backwards-compatibility to get
     * element access to a ghosted TrilinosWrappers::MPI::Vector inside the
     * library.
     */
    void
    reinit(const TrilinosWrappers::MPI::Vector &trilinos_vec);
#endif

    /**
     * Apply the functor @p func to each element of the vector. The functor
     * should look like
     * @code
     * struct Functor
     * {
     *   void operator() (Number &value);
     * };
     * @endcode
     *
     * @note This function requires that the header read_write_vector.templates.h
     * be included.
     */
    template <typename Functor>
    void
    apply(const Functor &func);

    /**
     * Swap the contents of this vector and the other vector @p v. One could
     * do this operation with a temporary variable and copying over the data
     * elements, but this function is significantly more efficient since it
     * only swaps the pointers to the data of the two vectors and therefore
     * does not need to allocate temporary storage and move data around.
     *
     * This function is analogous to the @p swap function of all C++
     * standard containers. Also, there is a global function
     * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
     * analogy to standard functions.
     */
    void
    swap(ReadWriteVector<Number> &v) noexcept;

    /**
     * Copies the data and the IndexSet of the input vector @p in_vector.
     */
    ReadWriteVector<Number> &
    operator=(const ReadWriteVector<Number> &in_vector);

    /**
     * Copies the data and the IndexSet of the input vector @p in_vector.
     */
    template <typename Number2>
    ReadWriteVector<Number> &
    operator=(const ReadWriteVector<Number2> &in_vector);

    /**
     * Sets all elements of the vector to the scalar @p s. This operation is
     * only allowed if @p s is equal to zero.
     */
    ReadWriteVector<Number> &
    operator=(const Number s);

    /**
     * Imports all the elements present in the vector's IndexSet from the
     * input vector @p vec. VectorOperation::values @p operation
     * is used to decide if the elements in @p V should be added to the
     * current vector or replace the current elements.
     *
     * @note The parameter @p communication_pattern is ignored since we are
     *   dealing with a serial vector here.
     */
    void
    import_elements(
      const dealii::Vector<Number> &vec,
      VectorOperation::values       operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    DEAL_II_DEPRECATED
    void
    import(const dealii::Vector<Number> &V,
           VectorOperation::values       operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }

    /**
     * Imports all the elements present in the vector's IndexSet from the
     * input vector @p vec. VectorOperation::values @p operation
     * is used to decide if the elements in @p V should be added to the
     * current vector or replace the current elements. The last parameter can
     * be used if the same communication pattern is used multiple times. This
     * can be used to improve performance.
     */
    template <typename MemorySpace>
    void
    import_elements(
      const distributed::Vector<Number, MemorySpace> &vec,
      VectorOperation::values                         operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    template <typename MemorySpace>
    DEAL_II_DEPRECATED void
    import(const distributed::Vector<Number, MemorySpace> &V,
           VectorOperation::values                         operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }


#ifdef DEAL_II_WITH_PETSC
    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p petsc_vec. VectorOperation::values @p operation is used to decide
     * if the elements in @p V should be added to the current vector or replace
     * the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    void
    import_elements(
      const PETScWrappers::MPI::Vector &petsc_vec,
      VectorOperation::values           operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    DEAL_II_DEPRECATED
    void
    import(const PETScWrappers::MPI::Vector &V,
           VectorOperation::values           operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }
#endif

#ifdef DEAL_II_WITH_TRILINOS
    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p trilinos_vec. VectorOperation::values @p operation is used to
     * decide if the elements in @p V should be added to the current vector or
     * replace the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     *
     * @note The @p trilinos_vec is not allowed to have ghost entries.
     */
    void
    import_elements(
      const TrilinosWrappers::MPI::Vector &trilinos_vec,
      VectorOperation::values              operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    DEAL_II_DEPRECATED
    void
    import(const TrilinosWrappers::MPI::Vector &V,
           VectorOperation::values              operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p tpetra_vec. VectorOperation::values @p operation is used to
     * decide if the elements in @p V should be added to the current vector or
     * replace the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    template <typename MemorySpace, typename Dummy = Number>
    std::enable_if_t<std::is_same_v<Dummy, Number> &&
                     dealii::is_tpetra_type<Number>::value>
    import_elements(
      const TpetraWrappers::Vector<Number, MemorySpace> &tpetra_vec,
      VectorOperation::values                            operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    template <typename MemorySpace, typename Dummy = Number>
    DEAL_II_DEPRECATED std::enable_if_t<std::is_same_v<Dummy, Number> &&
                                        dealii::is_tpetra_type<Number>::value>
    import(const TpetraWrappers::Vector<Number, MemorySpace> &V,
           VectorOperation::values                            operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }
#  endif

    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p epetra_vec. VectorOperation::values @p operation is used to
     * decide if the elements in @p V should be added to the current vector or
     * replace the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    void
    import_elements(
      const EpetraWrappers::Vector &epetra_vec,
      VectorOperation::values       operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern = {});

    /**
     * @deprecated Use import_elements() instead.
     */
    DEAL_II_DEPRECATED
    void
    import(const EpetraWrappers::Vector &V,
           VectorOperation::values       operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {})
    {
      import_elements(V, operation, communication_pattern);
    }
#endif

    /**
     * The value returned by this function denotes the dimension of the vector
     * spaces that are modeled by objects of this kind. However, objects of
     * the current class do not actually stores all elements of vectors of
     * this space but may, in fact store only a subset. The number of elements
     * stored is returned by n_elements() and is smaller or equal to the
     * number returned by the current function.
     */
    size_type
    size() const override;

    /**
     * Return the local size of the vector, i.e., the number of indices
     * owned locally.
     */
    size_type
    locally_owned_size() const;

    /**
     * Return the IndexSet that represents the indices of the elements stored.
     */
    const IndexSet &
    get_stored_elements() const;

    /**
     * Make the @p ReadWriteVector class a bit like the <tt>vector<></tt>
     * class of the C++ standard library by returning iterators to the start
     * and end of the <i>locally stored</i> elements of this vector.
     */
    iterator
    begin();

    /**
     * Return constant iterator to the start of the locally stored elements
     * of the vector.
     */
    const_iterator
    begin() const;

    /**
     * Return an iterator pointing to the element past the end of the array
     * of locally stored entries.
     */
    iterator
    end();

    /**
     * Return a constant iterator pointing to the element past the end of the
     * array of the locally stored entries.
     */
    const_iterator
    end() const;
    /** @} */


    /**
     * @name 2: Data-Access
     */
    /** @{ */

    /**
     * Read access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     */
    Number
    operator()(const size_type global_index) const;

    /**
     * Read and write access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     */
    Number &
    operator()(const size_type global_index);

    /**
     * Read access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     *
     * This function does the same thing as operator().
     */
    Number
    operator[](const size_type global_index) const;

    /**
     * Read and write access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     *
     * This function does the same thing as operator().
     */
    Number &
    operator[](const size_type global_index);

    /**
     * Instead of getting individual elements of a vector via operator(),
     * this function allows getting a whole set of elements at once. The
     * indices of the elements to be read are stated in the first argument, the
     * corresponding values are returned in the second.
     *
     * If the current vector is called @p v, then this function is the equivalent
     * to the code
     * @code
     *   for (unsigned int i=0; i<indices.size(); ++i)
     *     values[i] = v[indices[i]];
     * @endcode
     *
     * @pre The sizes of the @p indices and @p values arrays must be identical.
     */
    template <typename Number2>
    void
    extract_subvector_to(const std::vector<size_type> &indices,
                         std::vector<Number2>         &values) const;

    /**
     * Extract a range of elements all at once.
     */
    virtual void
    extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      const ArrayView<Number>                        &entries) const override;

    /**
     * Instead of getting individual elements of a vector via operator(),
     * this function allows getting a whole set of elements at once. In
     * contrast to the previous function, this function obtains the
     * indices of the elements by dereferencing all elements of the iterator
     * range provided by the first two arguments, and puts the vector
     * values into memory locations obtained by dereferencing a range
     * of iterators starting at the location pointed to by the third
     * argument.
     *
     * If the current vector is called @p v, then this function is the equivalent
     * to the code
     * @code
     *   ForwardIterator indices_p = indices_begin;
     *   OutputIterator  values_p  = values_begin;
     *   while (indices_p != indices_end)
     *   {
     *     *values_p = v[*indices_p];
     *     ++indices_p;
     *     ++values_p;
     *   }
     * @endcode
     *
     * @pre It must be possible to write into as many memory locations
     *   starting at @p values_begin as there are iterators between
     *   @p indices_begin and @p indices_end.
     */
    template <typename ForwardIterator, typename OutputIterator>
    void
    extract_subvector_to(ForwardIterator       indices_begin,
                         const ForwardIterator indices_end,
                         OutputIterator        values_begin) const;

    /**
     * Read access to the data field specified by @p local_index. When you
     * access elements in the order in which they are stored, it is necessary
     * that you know in which they are stored. In other words, you need to
     * know the map between the global indices of the elements this class
     * stores, and the local indices into the contiguous array of these global
     * elements. For this, see the general documentation of this class.
     *
     * Performance: Direct array access (fast).
     */
    Number
    local_element(const size_type local_index) const;

    /**
     * Read and write access to the data field specified by @p local_index.
     * When you access elements in the order in which they are stored, it is
     * necessary that you know in which they are stored. In other words, you
     * need to know the map between the global indices of the elements this
     * class stores, and the local indices into the contiguous array of these
     * global elements. For this, see the general documentation of this class.
     *
     * Performance: Direct array access (fast).
     */
    Number &
    local_element(const size_type local_index);
    /** @} */


    /**
     * @name 3: Modification of vectors
     */
    /** @{ */

    /**
     * This function adds a whole set of values stored in @p values to the
     * vector components specified by @p indices.
     */
    template <typename Number2>
    void
    add(const std::vector<size_type> &indices,
        const std::vector<Number2>   &values);

    /**
     * This function is similar to the previous one but takes a
     * ReadWriteVector of values.
     */
    template <typename Number2>
    void
    add(const std::vector<size_type>   &indices,
        const ReadWriteVector<Number2> &values);

    /**
     * Take an address where <tt>n_elements</tt> are stored contiguously and
     * add them into the vector. Handles all cases which are not covered by
     * the other two <tt>add()</tt> functions above.
     */
    template <typename Number2>
    void
    add(const size_type  n_elements,
        const size_type *indices,
        const Number2   *values);

    /**
     * Prints the vector to the output stream @p out.
     */
    void
    print(std::ostream      &out,
          const unsigned int precision  = 3,
          const bool         scientific = true) const;

    /**
     * Return the memory consumption of this class in bytes.
     */
    std::size_t
    memory_consumption() const;
    /** @} */

  protected:
#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * Import all the elements present in the vector's IndexSet from the input
     * vector @p tpetra_vector. This is an helper function and it should not be
     * used directly.
     */
    template <typename MemorySpace, typename Dummy = Number>
    std::enable_if_t<std::is_same_v<Dummy, Number> &&
                     dealii::is_tpetra_type<Number>::value>
    import_elements(
      const Tpetra::
        Vector<Number, int, types::signed_global_dof_index, MemorySpace>
                             &tpetra_vector,
      const IndexSet         &locally_owned_elements,
      VectorOperation::values operation,
      const MPI_Comm          mpi_comm,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern);
#  endif

    /**
     * Import all the elements present in the vector's IndexSet from the input
     * vector @p multivector. This is an helper function and it should not be
     * used directly.
     */
    void
    import_elements(
      const Epetra_MultiVector &multivector,
      const IndexSet           &locally_owned_elements,
      VectorOperation::values   operation,
      const MPI_Comm            mpi_comm,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern);
#endif

    /**
     * Return the local position of @p global_index.
     */
    unsigned int
    global_to_local(const types::global_dof_index global_index) const;

#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * Return a TpetraWrappers::CommunicationPattern and store it for future
     * use.
     */
    template <typename MemorySpace = dealii::MemorySpace::Host>
    TpetraWrappers::CommunicationPattern<MemorySpace>
    create_tpetra_comm_pattern(const IndexSet &source_index_set,
                               const MPI_Comm  mpi_comm);
#  endif

    /**
     * Return a EpetraWrappers::CommunicationPattern and store it for future
     * use.
     */
    EpetraWrappers::CommunicationPattern
    create_epetra_comm_pattern(const IndexSet &source_index_set,
                               const MPI_Comm  mpi_comm);
#endif

    /**
     * Indices of the elements stored.
     */
    IndexSet stored_elements;

    /**
     * IndexSet of the elements of the last imported vector;
     */
    IndexSet source_stored_elements;

    /**
     * CommunicationPattern for the communication between the
     * source_stored_elements IndexSet and the current vector.
     */
    std::shared_ptr<Utilities::MPI::CommunicationPatternBase> comm_pattern;

    /**
     * Locally stored elements.
     */
    AlignedVector<Number> values;

    /**
     * For parallel loops with TBB, this member variable stores the affinity
     * information of loops.
     */
    mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
      thread_loop_partitioner;

    // Make all other ReadWriteVector types friends.
    template <typename Number2>
    friend class ReadWriteVector;

  private:
    /**
     * This class provides a wrapper around a Functor which acts on
     * single elements of the vector. This is necessary to use
     * tbb::parallel_for which requires a TBBForFunctor.
     */
    template <typename Functor>
    class FunctorTemplate
    {
    public:
      /**
       * Constructor. Take a functor and store a copy of it.
       */
      FunctorTemplate(ReadWriteVector<Number> &parent, const Functor &functor);

      /**
       * Evaluate the element with the stored copy of the functor.
       */
      virtual void
      operator()(const size_type begin, const size_type end);

    private:
      /**
       * Alias to the ReadWriteVector object that owns the FunctorTemplate.
       */
      ReadWriteVector &parent;

      /**
       * Copy of the functor.
       */
      const Functor &functor;
    };
  };

  /** @} */


  /*---------------------------- Inline functions ---------------------------*/

#ifndef DOXYGEN


  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector()
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(0, true);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(
    const ReadWriteVector<Number> &v)
  {
    this->operator=(v);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(const size_type size)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(size, false);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(
    const IndexSet &locally_stored_indices)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(locally_stored_indices);
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::size() const
  {
    return stored_elements.size();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::locally_owned_size() const
  {
    return stored_elements.n_elements();
  }



  template <typename Number>
  inline const IndexSet &
  ReadWriteVector<Number>::get_stored_elements() const
  {
    return stored_elements;
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::begin()
  {
    return values.begin();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::begin() const
  {
    return values.begin();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::end()
  {
    return values.end();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::end() const
  {
    return values.end();
  }



  template <typename Number>
  inline Number
  ReadWriteVector<Number>::operator()(const size_type global_index) const
  {
    return values[global_to_local(global_index)];
  }



  template <typename Number>
  inline Number &
  ReadWriteVector<Number>::operator()(const size_type global_index)
  {
    return values[global_to_local(global_index)];
  }



  template <typename Number>
  inline Number
  ReadWriteVector<Number>::operator[](const size_type global_index) const
  {
    return operator()(global_index);
  }



  template <typename Number>
  inline Number &
  ReadWriteVector<Number>::operator[](const size_type global_index)
  {
    return operator()(global_index);
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::extract_subvector_to(
    const std::vector<size_type> &indices,
    std::vector<Number2>         &extracted_values) const
  {
    for (size_type i = 0; i < indices.size(); ++i)
      extracted_values[i] = operator()(indices[i]);
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::extract_subvector_to(
    const ArrayView<const types::global_dof_index> &indices,
    const ArrayView<Number>                        &entries) const
  {
    AssertDimension(indices.size(), entries.size());
    for (unsigned int i = 0; i < indices.size(); ++i)
      entries[i] = (*this)[indices[i]];
  }



  template <typename Number>
  template <typename ForwardIterator, typename OutputIterator>
  inline void
  ReadWriteVector<Number>::extract_subvector_to(
    ForwardIterator       indices_begin,
    const ForwardIterator indices_end,
    OutputIterator        values_begin) const
  {
    while (indices_begin != indices_end)
      {
        *values_begin = operator()(*indices_begin);
        ++indices_begin;
        ++values_begin;
      }
  }



  template <typename Number>
  inline Number
  ReadWriteVector<Number>::local_element(const size_type local_index) const
  {
    AssertIndexRange(local_index, this->locally_owned_size());

    return values[local_index];
  }



  template <typename Number>
  inline Number &
  ReadWriteVector<Number>::local_element(const size_type local_index)
  {
    AssertIndexRange(local_index, this->locally_owned_size());

    return values[local_index];
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const std::vector<size_type> &indices,
                               const std::vector<Number2>   &values)
  {
    AssertDimension(indices.size(), values.size());
    add(indices.size(), indices.data(), values.data());
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const std::vector<size_type>   &indices,
                               const ReadWriteVector<Number2> &values)
  {
    const size_type size = indices.size();
    for (size_type i = 0; i < size; ++i)
      {
        Assert(
          numbers::is_finite(values[i]),
          ExcMessage(
            "The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values[indices[i]];
      }
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const size_type  n_indices,
                               const size_type *indices,
                               const Number2   *values_to_add)
  {
    for (size_type i = 0; i < n_indices; ++i)
      {
        Assert(
          numbers::is_finite(values[i]),
          ExcMessage(
            "The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values_to_add[i];
      }
  }



  template <typename Number>
  inline unsigned int
  ReadWriteVector<Number>::global_to_local(
    const types::global_dof_index global_index) const
  {
    // the following will throw an exception if the global_index is not
    // in the remaining_elements
    return static_cast<unsigned int>(
      stored_elements.index_within_set(global_index));
  }



  template <typename Number>
  template <typename Functor>
  inline ReadWriteVector<Number>::FunctorTemplate<Functor>::FunctorTemplate(
    ReadWriteVector<Number> &parent,
    const Functor           &functor)
    : parent(parent)
    , functor(functor)
  {}



  template <typename Number>
  template <typename Functor>
  void
  ReadWriteVector<Number>::FunctorTemplate<Functor>::operator()(
    const size_type begin,
    const size_type end)
  {
    for (size_type i = begin; i < end; ++i)
      functor(parent.values[i]);
  }

#endif // ifndef DOXYGEN

} // end of namespace LinearAlgebra



/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relatesalso Vector
 */
template <typename Number>
inline void
swap(LinearAlgebra::ReadWriteVector<Number> &u,
     LinearAlgebra::ReadWriteVector<Number> &v) noexcept
{
  u.swap(v);
}


DEAL_II_NAMESPACE_CLOSE

#endif
