// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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


#ifndef dealii_aligned_vector_h
#define dealii_aligned_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>

// boost::serialization::make_array used to be in array.hpp, but was
// moved to a different file in BOOST 1.64
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/serialization/split_member.hpp>

#include <cstring>
#include <memory>
#include <type_traits>



DEAL_II_NAMESPACE_OPEN


/**
 * This is a replacement class for std::vector to be used in combination with
 * VectorizedArray and derived data types. It allocates memory aligned to
 * addresses of a vectorized data type (in order to avoid segmentation faults
 * when a variable of type VectorizedArray which the compiler assumes to be
 * aligned to certain memory addresses does not actually follow these rules).
 * This could also be achieved by proving std::vector with a user-defined
 * allocator. On the other hand, writing an own small vector class lets us
 * implement parallel copy and move operations with TBB, insert deal.II-style
 * assertions, and cut some unnecessary functionality. Note that this vector
 * is a bit more memory-consuming than std::vector because of alignment, so it
 * is recommended to only use this vector on long vectors.
 *
 * @p author Katharina Kormann, Martin Kronbichler, 2011
 */
template <class T>
class AlignedVector
{
public:
  /**
   * Declare standard types used in all containers. These types parallel those
   * in the <tt>C++</tt> standard libraries <tt>vector<...></tt> class.
   */
  using value_type      = T;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;
  using size_type       = std::size_t;

  /**
   * Empty constructor. Sets the vector size to zero.
   */
  AlignedVector();

  /**
   * Set the vector size to the given size and initializes all elements with
   * T().
   *
   * @dealiiOperationIsMultithreaded
   */
  explicit AlignedVector(const size_type size, const T &init = T());

  /**
   * Destructor.
   */
  ~AlignedVector();

  /**
   * Copy constructor.
   *
   * @dealiiOperationIsMultithreaded
   */
  AlignedVector(const AlignedVector<T> &vec);

  /**
   * Move constructor. Create a new aligned vector by stealing the contents of
   * @p vec.
   */
  AlignedVector(AlignedVector<T> &&vec) noexcept;

  /**
   * Assignment to the input vector @p vec.
   *
   * @dealiiOperationIsMultithreaded
   */
  AlignedVector &
  operator=(const AlignedVector<T> &vec);

  /**
   * Move assignment operator.
   */
  AlignedVector &
  operator=(AlignedVector<T> &&vec) noexcept;

  /**
   * Change the size of the vector. If the new size is larger than the
   * previous size, then new elements will be added to the end of the
   * vector; these elements will remain uninitialized (i.e., left in
   * an undefined state) if `std::is_trivial<T>` is `true`, and
   * will be default initialized if `std::is_trivial<T>` is `false`.
   * See [here](https://en.cppreference.com/w/cpp/types/is_trivial) for
   * a definition of what `std::is_trivial` does.
   *
   * If the new size is less than the previous size, then the last few
   * elements will be destroyed if `std::is_trivial<T>` will be `false`
   * or will simply be ignored in the future if
   * `std::is_trivial<T>` is `true`.
   *
   * As a consequence of the outline above, the "_fast" suffix of this
   * function refers to the fact that for "trivial" classes `T`, this
   * function omits constructor/destructor calls and in particular the
   * initialization of new elements.
   *
   * @note This method can only be invoked for classes @p T that define a
   * default constructor, @p T(). Otherwise, compilation will fail.
   */
  void
  resize_fast(const size_type new_size);

  /**
   * Change the size of the vector. It keeps old elements previously
   * available, and initializes each newly added element to a
   * default-constructed object of type @p T.
   *
   * If the new vector size is shorter than the old one, the memory is
   * not immediately released unless the new size is zero; however,
   * the size of the current object is of course set to the requested
   * value. The destructors of released elements are also called.
   *
   * @dealiiOperationIsMultithreaded
   */
  void
  resize(const size_type new_size);

  /**
   * Change the size of the vector. It keeps old elements previously
   * available, and initializes each newly added element with the
   * provided initializer.
   *
   * If the new vector size is shorter than the old one, the memory is
   * not immediately released unless the new size is zero; however,
   * the size of the current object is of course set to the requested
   * value.
   *
   * @note This method can only be invoked for classes that define the copy
   * assignment operator. Otherwise, compilation will fail.
   *
   * @dealiiOperationIsMultithreaded
   */
  void
  resize(const size_type new_size, const T &init);

  /**
   * Reserve memory space for @p new_allocated_size elements.
   *
   * If the argument @p new_allocated_size is less than the current number of stored
   * elements (as indicated by calling size()), then this function does not
   * do anything at all. Except if the argument @p new_allocated_size is set
   * to zero, then all previously allocated memory is released (this operation
   * then being equivalent to directly calling the clear() function).
   *
   * In order to avoid too frequent reallocation (which involves copy of the
   * data), this function doubles the amount of memory occupied when the given
   * size is larger than the previously allocated size.
   *
   * Note that this function only changes the amount of elements the object
   * *can* store, but not the number of elements it *actually* stores. As
   * a consequence, no constructors or destructors of newly created objects
   * are run, though the existing elements may be moved to a new location (which
   * involves running the move constructor at the new location and the
   * destructor at the old location).
   */
  void
  reserve(const size_type new_allocated_size);

  /**
   * Releases all previously allocated memory and leaves the vector in a state
   * equivalent to the state after the default constructor has been called.
   */
  void
  clear();

  /**
   * Inserts an element at the end of the vector, increasing the vector size
   * by one. Note that the allocated size will double whenever the previous
   * space is not enough to hold the new element.
   */
  void
  push_back(const T in_data);

  /**
   * Return the last element of the vector (read and write access).
   */
  reference
  back();

  /**
   * Return the last element of the vector (read-only access).
   */
  const_reference
  back() const;

  /**
   * Inserts several elements at the end of the vector given by a range of
   * elements.
   */
  template <typename ForwardIterator>
  void
  insert_back(ForwardIterator begin, ForwardIterator end);

  /**
   * Fills the vector with size() copies of a default constructed object.
   *
   * @note Unlike the other fill() function, this method can also be
   * invoked for classes that do not define a copy assignment
   * operator.
   *
   * @dealiiOperationIsMultithreaded
   */
  void
  fill();

  /**
   * Fills the vector with size() copies of the given input.
   *
   * @note This method can only be invoked for classes that define the copy
   * assignment operator. Otherwise, compilation will fail.
   *
   * @dealiiOperationIsMultithreaded
   */
  void
  fill(const T &element);

  /**
   * This function replicates the state found on the process indicated by
   * @p root_process across all processes of the MPI communicator. The current
   * state found on any of the processes other than @p root_process is lost
   * in this process. One can imagine this operation to act like a call to
   * Utilities::MPI::broadcast() from the root process to all other processes,
   * though in practice the function may try to move the data into shared
   * memory regions on each of the machines that host MPI processes and
   * let all MPI processes on this machine then access this shared memory
   * region instead of keeping their own copy.
   *
   * The intent of this function is to quickly exchange large arrays from
   * one process to others, rather than having to compute or create it on
   * all processes. This is specifically the case for data loaded from
   * disk -- say, large data tables -- that are more easily dealt with by
   * reading once and then distributing across all processes in an MPI
   * universe, than letting each process read the data from disk itself.
   * Specifically, the use of shared memory regions allows for replicating
   * the data only once per multicore machine in the MPI universe, rather
   * than replicating data once for each MPI process. This results in
   * large memory savings if the data is large on today's machines that
   * can easily house several dozen MPI processes per shared memory
   * space. This use case is outlined in the TableBase class documentation
   * as the current function is called from
   * TableBase::replicate_across_communicator(). Indeed, the primary rationale
   * for this function is to enable sharing data tables based on TableBase
   * across MPI processes.
   *
   * This function does not imply a model of keeping data on different processes
   * in sync, as parallel::distributed::Vector and other vector classes do where
   * there exists a notion of certain elements of the vector owned by each
   * process and possibly ghost elements that are mirrored from its owning
   * process to other processes. Rather, the elements of the current object are
   * simply copied to the other processes, and it is useful to think of this
   * operation as creating a set of `const` AlignedVector objects on all
   * processes that should not be changed any more after the replication
   * operation, as this is the only way to ensure that the vectors remain the
   * same on all processes. This is particularly true because of the use of
   * shared memory regions where any modification of a vector element on one MPI
   * process may also result in a modification of elements visible on other
   * processes, assuming they are located within one shared memory node.
   *
   * @note The use of shared memory between MPI processes requires
   *   that the detected MPI installation supports the necessary operations.
   *   This is the case for MPI 3.0 and higher.
   *
   * @note This function is not cheap. It needs to create sub-communicators
   *   of the provided @p communicator object, which is generally an expensive
   *   operation. Likewise, the generation of shared memory spaces is not
   *   a cheap operation. As a consequence, this function primarily makes
   *   sense when the goal is to share large read-only data tables among
   *   processes; examples are data tables that are loaded at start-up
   *   time and then used over the course of the run time of the program.
   *   In such cases, the start-up cost of running this function can be
   *   amortized over time, and the potential memory savings from not having to
   *   store the table on each process may be substantial on machines with
   *   large core counts on which many MPI processes run on the same machine.
   *
   * @note This function only makes sense if the data type `T` is
   *   "self-contained", i.e., all if its information is stored in its
   *   member variables, and if none of the member variables are pointers
   *   to other parts of the memory. This is because if a type `T` does
   *   have pointers to other parts of memory, then moving `T` into
   *   a shared memory space does not result in the other processes having
   *   access to data that the object points to with its member variable
   *   pointers: These continue to live only on one process, and are
   *   typically in memory areas not accessible to the other processes.
   *   As a consequence, the usual use case for this function is to share
   *   arrays of simple objects such as `double`s or `int`s.
   *
   * @note After calling this function, objects on different MPI processes
   *   share a common state. That means that certain operations become
   *   "collective", i.e., they must be called on all participating
   *   processors at the same time. In particular, you can no longer call
   *   resize(), reserve(), or clear() on one MPI process -- you have to do
   *   so on all processes at the same time, because they have to communicate
   *   for these operations. If you do not do so, you will likely get
   *   a deadlock that may be difficult to debug. By extension, this rule of
   *   only collectively resizing extends to this function itself: You can
   *   not call it twice in a row because that implies that first all but the
   *   `root_process` throw away their data, which is not a collective
   *   operation. Generally, these restrictions on what can and can not be
   *   done hint at the correctness of the comments above: You should treat
   *   an AlignedVector on which the current function has been called as
   *   `const`, on which no further operations can be performed until
   *   the destructor is called.
   */
  void
  replicate_across_communicator(const MPI_Comm &   communicator,
                                const unsigned int root_process);

  /**
   * Swaps the given vector with the calling vector.
   */
  void
  swap(AlignedVector<T> &vec);

  /**
   * Return whether the vector is empty, i.e., its size is zero.
   */
  bool
  empty() const;

  /**
   * Return the size of the vector.
   */
  size_type
  size() const;

  /**
   * Return the capacity of the vector, i.e., the size this vector can hold
   * without reallocation. Note that capacity() >= size().
   */
  size_type
  capacity() const;

  /**
   * Read-write access to entry @p index in the vector.
   */
  reference operator[](const size_type index);

  /**
   * Read-only access to entry @p index in the vector.
   */
  const_reference operator[](const size_type index) const;

  /**
   * Return a pointer to the underlying data buffer.
   */
  pointer
  data();

  /**
   * Return a const pointer to the underlying data buffer.
   */
  const_pointer
  data() const;

  /**
   * Return a read and write pointer to the beginning of the data array.
   */
  iterator
  begin();

  /**
   * Return a read and write pointer to the end of the data array.
   */
  iterator
  end();

  /**
   * Return a read-only pointer to the beginning of the data array.
   */
  const_iterator
  begin() const;

  /**
   * Return a read-only pointer to the end of the data array.
   */
  const_iterator
  end() const;

  /**
   * Return the memory consumption of the allocated memory in this class. If
   * the underlying type @p T allocates memory by itself, this memory is not
   * counted.
   */
  size_type
  memory_consumption() const;

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * Write and read the data of this object from a stream for the purpose
   * of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

private:
  /**
   * Pointer to actual data array.
   */
  std::unique_ptr<T[], std::function<void(T *)>> elements;

  /**
   * Pointer to one past the last valid value.
   */
  T *used_elements_end;

  /**
   * Pointer to the end of the allocated memory.
   */
  T *allocated_elements_end;
};


// ------------------------------- inline functions --------------------------

/**
 * This namespace defines the copy and set functions used in AlignedVector.
 * These functions operate in parallel when there are enough elements in the
 * vector.
 */
namespace internal
{
  /**
   * A class that actually issues the copy commands in AlignedVector.
   * This class is based on the specialized for loop base class
   * ParallelForLoop in parallel.h whose purpose is the following: When
   * calling a parallel for loop on AlignedVector with apply_to_subranges, it
   * generates different code for every different argument we might choose (as
   * it is templated). This gives a lot of code (e.g. it triples the memory
   * required for compiling the file matrix_free.cc and the final object size
   * is several times larger) which is completely useless. Therefore, this
   * class channels all copy commands through one call to apply_to_subrange
   * for all possible types, which makes the copy operation much cleaner
   * (thanks to a virtual function, whose cost is negligible in this context).
   *
   * @relatesalso AlignedVector
   */
  template <typename T>
  class AlignedVectorCopy : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise works in serial. Copies the data from the half-open
     * interval between @p source_begin and @p source_end to array starting at
     * @p destination (by calling the copy constructor with placement new).
     *
     * The elements from the source array are simply copied via the placement
     * new copy constructor.
     */
    AlignedVectorCopy(const T *const source_begin,
                      const T *const source_end,
                      T *const       destination)
      : source_(source_begin)
      , destination_(destination)
    {
      Assert(source_end >= source_begin, ExcInternalError());
      Assert(source_end == source_begin || destination != nullptr,
             ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        AlignedVectorCopy::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * This method moves elements from the source to the destination given in
     * the constructor on a subrange given by two integers.
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      if (end == begin)
        return;

      // for classes trivial assignment can use memcpy. cast element to
      // (void*) to silence compiler warning for virtual classes (they will
      // never arrive here because they are non-trivial).

      if (std::is_trivial<T>::value == true)
        std::memcpy(static_cast<void *>(destination_ + begin),
                    static_cast<const void *>(source_ + begin),
                    (end - begin) * sizeof(T));
      else
        for (std::size_t i = begin; i < end; ++i)
          new (&destination_[i]) T(source_[i]);
    }

  private:
    const T *const source_;
    T *const       destination_;
  };


  /**
   * Like AlignedVectorCopy, but use move operations instead of copy operations.
   *
   * @relatesalso AlignedVector
   */
  template <typename T>
  class AlignedVectorMove : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise works in serial. Moves the data from the half-open
     * interval between @p source_begin and @p source_end to array starting at
     * @p destination (by calling the move constructor with placement new).
     *
     * The data is moved between the two arrays by invoking the destructor on
     * the source range (preparing for a subsequent call to free).
     */
    AlignedVectorMove(T *const source_begin,
                      T *const source_end,
                      T *const destination)
      : source_(source_begin)
      , destination_(destination)
    {
      Assert(source_end >= source_begin, ExcInternalError());
      Assert(source_end == source_begin || destination != nullptr,
             ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        AlignedVectorMove::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * This method moves elements from the source to the destination given in
     * the constructor on a subrange given by two integers.
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      if (end == begin)
        return;

      // for classes trivial assignment can use memcpy. cast element to
      // (void*) to silence compiler warning for virtual classes (they will
      // never arrive here because they are non-trivial).

      if (std::is_trivial<T>::value == true)
        std::memcpy(static_cast<void *>(destination_ + begin),
                    static_cast<void *>(source_ + begin),
                    (end - begin) * sizeof(T));
      else
        for (std::size_t i = begin; i < end; ++i)
          {
            // initialize memory (copy construct by placement new), and
            // destruct the source
            new (&destination_[i]) T(std::move(source_[i]));
            source_[i].~T();
          }
    }

  private:
    T *const source_;
    T *const destination_;
  };


  /**
   * Class that issues the set commands for AlignedVector.
   *
   * @tparam initialize_memory Sets whether the set command should
   * initialize memory (with a call to the copy constructor) or rather use the
   * copy assignment operator. A template is necessary to select the
   * appropriate operation since some classes might define only one of those
   * two operations.
   *
   * @relatesalso AlignedVector
   */
  template <typename T, bool initialize_memory>
  class AlignedVectorSet : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise work in serial.
     */
    AlignedVectorSet(const std::size_t size,
                     const T &         element,
                     T *const          destination)
      : element_(element)
      , destination_(destination)
      , trivial_element(false)
    {
      if (size == 0)
        return;
      Assert(destination != nullptr, ExcInternalError());

      // do not use memcmp for long double because on some systems it does not
      // completely fill its memory and may lead to false positives in
      // e.g. valgrind
      if (std::is_trivial<T>::value == true &&
          std::is_same<T, long double>::value == false)
        {
          const unsigned char zero[sizeof(T)] = {};
          // cast element to (void*) to silence compiler warning for virtual
          // classes (they will never arrive here because they are
          // non-trivial).
          if (std::memcmp(zero,
                          static_cast<const void *>(&element),
                          sizeof(T)) == 0)
            trivial_element = true;
        }
      if (size < minimum_parallel_grain_size)
        AlignedVectorSet::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * This sets elements on a subrange given by two integers.
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      // for classes with trivial assignment of zero can use memset. cast
      // element to (void*) to silence compiler warning for virtual
      // classes (they will never arrive here because they are
      // non-trivial).
      if (std::is_trivial<T>::value == true && trivial_element)
        std::memset(static_cast<void *>(destination_ + begin),
                    0,
                    (end - begin) * sizeof(T));
      else
        copy_construct_or_assign(
          begin, end, std::integral_constant<bool, initialize_memory>());
    }

  private:
    const T &  element_;
    mutable T *destination_;
    bool       trivial_element;

    // copy assignment operation
    void
    copy_construct_or_assign(const std::size_t begin,
                             const std::size_t end,
                             std::integral_constant<bool, false>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        destination_[i] = element_;
    }

    // copy constructor (memory initialization)
    void
    copy_construct_or_assign(const std::size_t begin,
                             const std::size_t end,
                             std::integral_constant<bool, true>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        new (&destination_[i]) T(element_);
    }
  };



  /**
   * Class that issues the set commands for AlignedVector.
   *
   * @tparam initialize_memory Sets whether the set command should
   * initialize memory (with a call to the copy constructor) or rather use the
   * copy assignment operator. A template is necessary to select the
   * appropriate operation since some classes might define only one of those
   * two operations.
   *
   * @relatesalso AlignedVector
   */
  template <typename T, bool initialize_memory>
  class AlignedVectorDefaultInitialize
    : private dealii::parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size =
      160000 / sizeof(T) + 1;

  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise work in serial.
     */
    AlignedVectorDefaultInitialize(const std::size_t size, T *const destination)
      : destination_(destination)
    {
      if (size == 0)
        return;
      Assert(destination != nullptr, ExcInternalError());

      if (size < minimum_parallel_grain_size)
        AlignedVectorDefaultInitialize::apply_to_subrange(0, size);
      else
        apply_parallel(0, size, minimum_parallel_grain_size);
    }

    /**
     * This initializes elements on a subrange given by two integers.
     */
    virtual void
    apply_to_subrange(const std::size_t begin,
                      const std::size_t end) const override
    {
      // for classes with trivial assignment of zero can use memset. cast
      // element to (void*) to silence compiler warning for virtual
      // classes (they will never arrive here because they are
      // non-trivial).
      if (std::is_trivial<T>::value == true)
        std::memset(static_cast<void *>(destination_ + begin),
                    0,
                    (end - begin) * sizeof(T));
      else
        default_construct_or_assign(
          begin, end, std::integral_constant<bool, initialize_memory>());
    }

  private:
    mutable T *destination_;

    // copy assignment operation
    void
    default_construct_or_assign(const std::size_t begin,
                                const std::size_t end,
                                std::integral_constant<bool, false>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        destination_[i] = std::move(T());
    }

    // copy constructor (memory initialization)
    void
    default_construct_or_assign(const std::size_t begin,
                                const std::size_t end,
                                std::integral_constant<bool, true>) const
    {
      for (std::size_t i = begin; i < end; ++i)
        new (&destination_[i]) T;
    }
  };

} // end of namespace internal


#ifndef DOXYGEN


template <class T>
inline AlignedVector<T>::AlignedVector()
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{}



template <class T>
inline AlignedVector<T>::AlignedVector(const size_type size, const T &init)
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{
  if (size > 0)
    resize(size, init);
}



template <class T>
inline AlignedVector<T>::~AlignedVector()
{
  clear();
}



template <class T>
inline AlignedVector<T>::AlignedVector(const AlignedVector<T> &vec)
  : elements(nullptr, [](T *) { Assert(false, ExcInternalError()); })
  , used_elements_end(nullptr)
  , allocated_elements_end(nullptr)
{
  // copy the data from vec
  reserve(vec.size());
  used_elements_end = allocated_elements_end;
  internal::AlignedVectorCopy<T>(vec.elements.get(),
                                 vec.used_elements_end,
                                 elements.get());
}



template <class T>
inline AlignedVector<T>::AlignedVector(AlignedVector<T> &&vec) noexcept
  : elements(std::move(vec.elements))
  , used_elements_end(vec.used_elements_end)
  , allocated_elements_end(vec.allocated_elements_end)
{
  vec.elements               = nullptr;
  vec.used_elements_end      = nullptr;
  vec.allocated_elements_end = nullptr;
}



template <class T>
inline AlignedVector<T> &
AlignedVector<T>::operator=(const AlignedVector<T> &vec)
{
  resize(0);
  resize_fast(vec.used_elements_end - vec.elements.get());
  internal::AlignedVectorCopy<T>(vec.elements.get(),
                                 vec.used_elements_end,
                                 elements.get());
  return *this;
}



template <class T>
inline AlignedVector<T> &
AlignedVector<T>::operator=(AlignedVector<T> &&vec) noexcept
{
  clear();

  // Move the actual data
  elements = std::move(vec.elements);

  // Then also steal the other pointers and clear them in the original object:
  used_elements_end      = vec.used_elements_end;
  allocated_elements_end = vec.allocated_elements_end;

  vec.used_elements_end      = nullptr;
  vec.allocated_elements_end = nullptr;

  return *this;
}



template <class T>
inline void
AlignedVector<T>::resize_fast(const size_type new_size)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released. doing it backward
      // releases the elements in reverse order as compared to how they were
      // created
      if (std::is_trivial<T>::value == false)
        {
          while (used_elements_end != elements.get() + new_size)
            (--used_elements_end)->~T();
        }
      else
        used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // need to still set the values in case the class is non-trivial because
      // virtual classes etc. need to run their (default) constructor
      if (std::is_trivial<T>::value == false)
        dealii::internal::AlignedVectorDefaultInitialize<T, true>(
          new_size - old_size, elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::resize(const size_type new_size)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released. doing it backward
      // releases the elements in reverse order as compared to how they were
      // created
      if (std::is_trivial<T>::value == false)
        {
          while (used_elements_end != elements.get() + new_size)
            (--used_elements_end)->~T();
        }
      else
        used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // finally set the values to the default initializer
      dealii::internal::AlignedVectorDefaultInitialize<T, true>(
        new_size - old_size, elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::resize(const size_type new_size, const T &init)
{
  const size_type old_size = size();

  if (new_size == 0)
    clear();
  else if (new_size == old_size)
    {} // nothing to do here
  else if (new_size < old_size)
    {
      // call destructor on fields that are released. doing it backward
      // releases the elements in reverse order as compared to how they were
      // created
      if (std::is_trivial<T>::value == false)
        {
          while (used_elements_end != elements.get() + new_size)
            (--used_elements_end)->~T();
        }
      else
        used_elements_end = elements.get() + new_size;
    }
  else // new_size > old_size
    {
      // Allocate more space, and claim that space as used
      reserve(new_size);
      used_elements_end = elements.get() + new_size;

      // finally set the desired init values
      dealii::internal::AlignedVectorSet<T, true>(new_size - old_size,
                                                  init,
                                                  elements.get() + old_size);
    }
}



template <class T>
inline void
AlignedVector<T>::reserve(const size_type new_allocated_size)
{
  const size_type old_size           = used_elements_end - elements.get();
  const size_type old_allocated_size = allocated_elements_end - elements.get();
  if (new_allocated_size > old_allocated_size)
    {
      // if we continuously increase the size of the vector, we might be
      // reallocating a lot of times. therefore, try to increase the size more
      // aggressively
      const size_type new_size =
        std::max(new_allocated_size, 2 * old_allocated_size);

      // allocate and align along 64-byte boundaries (this is enough for all
      // levels of vectorization currently supported by deal.II)
      T *new_data_ptr;
      Utilities::System::posix_memalign(
        reinterpret_cast<void **>(&new_data_ptr), 64, new_size * sizeof(T));
      std::unique_ptr<T[], void (*)(T *)> new_data(new_data_ptr, [](T *ptr) {
        std::free(ptr);
      });

      // copy whatever elements we need to retain
      if (new_allocated_size > 0)
        dealii::internal::AlignedVectorMove<T>(elements.get(),
                                               elements.get() + old_size,
                                               new_data.get());

      // Now reset all of the member variables of the current object
      // based on the allocation above. Assigning to a std::unique_ptr
      // object also releases the previously pointed to memory.
      elements               = std::move(new_data);
      used_elements_end      = elements.get() + old_size;
      allocated_elements_end = elements.get() + new_size;
    }
  else if (new_allocated_size == 0)
    clear();
  else // size_alloc < allocated_size
    {} // nothing to do here
}



template <class T>
inline void
AlignedVector<T>::clear()
{
  if (elements != nullptr)
    {
      if (std::is_trivial<T>::value == false)
        while (used_elements_end != elements.get())
          (--used_elements_end)->~T();
    }
  elements               = nullptr;
  used_elements_end      = nullptr;
  allocated_elements_end = nullptr;
}



template <class T>
inline void
AlignedVector<T>::push_back(const T in_data)
{
  Assert(used_elements_end <= allocated_elements_end, ExcInternalError());
  if (used_elements_end == allocated_elements_end)
    reserve(std::max(2 * capacity(), static_cast<size_type>(16)));
  if (std::is_trivial<T>::value == false)
    new (used_elements_end++) T(in_data);
  else
    *used_elements_end++ = in_data;
}



template <class T>
inline typename AlignedVector<T>::reference
AlignedVector<T>::back()
{
  AssertIndexRange(0, size());
  T *field = used_elements_end - 1;
  return *field;
}



template <class T>
inline typename AlignedVector<T>::const_reference
AlignedVector<T>::back() const
{
  AssertIndexRange(0, size());
  const T *field = used_elements_end - 1;
  return *field;
}



template <class T>
template <typename ForwardIterator>
inline void
AlignedVector<T>::insert_back(ForwardIterator begin, ForwardIterator end)
{
  const unsigned int old_size = size();
  reserve(old_size + (end - begin));
  for (; begin != end; ++begin, ++used_elements_end)
    {
      if (std::is_trivial<T>::value == false)
        new (used_elements_end) T;
      *used_elements_end = *begin;
    }
}



template <class T>
inline void
AlignedVector<T>::fill()
{
  dealii::internal::AlignedVectorDefaultInitialize<T, false>(size(),
                                                             elements.get());
}



template <class T>
inline void
AlignedVector<T>::fill(const T &value)
{
  dealii::internal::AlignedVectorSet<T, false>(size(), value, elements.get());
}



template <class T>
inline void
AlignedVector<T>::replicate_across_communicator(const MPI_Comm &   communicator,
                                                const unsigned int root_process)
{
#  ifdef DEAL_II_WITH_MPI
#    if DEAL_II_MPI_VERSION_GTE(3, 0)

  // **** Step 0 ****
  // All but the root process no longer need their data, so release the memory
  // used to store the previous elements.
  if (Utilities::MPI::this_mpi_process(communicator) != root_process)
    {
      elements.reset();
      used_elements_end      = nullptr;
      allocated_elements_end = nullptr;
    }

  // **** Step 1 ****
  // Create communicators for each group of processes that can use
  // shared memory areas. Within each of these groups, we don't care about
  // which rank each of the old processes gets except that we would like to
  // make sure that the (global) root process will have rank=0 within
  // its own sub-communicator. We can do that through the third argument of
  // MPI_Comm_split_type (the "key") which is an integer meant to indicate the
  // order of processes within the split communicators, and we should set it to
  // zero for the root processes and one for all others -- which means that
  // for all of these other processes, MPI can choose whatever order it
  // wants because they have the same key (MPI then documents that these ties
  // will be broken according to these processes' rank in the old group).
  //
  // At least that's the theory. In practice, the MPI implementation where
  // this function was developed on does not seem to do that. (Bug report
  // is here: https://github.com/open-mpi/ompi/issues/8854)
  // We work around this by letting MPI_Comm_split_type choose whatever
  // rank it wants, and then reshuffle with MPI_Comm_split in a second
  // step -- not elegant, nor efficient, but seems to work:
  MPI_Comm shmem_group_communicator;
  {
    MPI_Comm shmem_group_communicator_temp;
    int      ierr = MPI_Comm_split_type(communicator,
                                   MPI_COMM_TYPE_SHARED,
                                   /* key */ 0,
                                   MPI_INFO_NULL,
                                   &shmem_group_communicator_temp);

    AssertThrowMPI(ierr);
    const int key =
      (Utilities::MPI::this_mpi_process(communicator) == root_process ? 0 : 1);
    ierr = MPI_Comm_split(shmem_group_communicator_temp,
                          /* color */ 0,
                          key,
                          &shmem_group_communicator);
    AssertThrowMPI(ierr);

    // Verify the explanation from above
    if (Utilities::MPI::this_mpi_process(communicator) == root_process)
      Assert(Utilities::MPI::this_mpi_process(shmem_group_communicator) == 0,
             ExcInternalError());

    // And get rid of the temporary communicator
    ierr = MPI_Comm_free(&shmem_group_communicator_temp);
  }
  const bool is_shmem_root =
    Utilities::MPI::this_mpi_process(shmem_group_communicator) == 0;

  // **** Step 2 ****
  // We then have to send the state of the current object from the
  // root process to one exemplar in each shmem group. To this end,
  // we create another subcommunicator that includes the ranks zero
  // of all shmem groups, and because of the trick above, we know
  // that this also includes the original root process.
  //
  // There are different ways of creating a "shmem_roots_communicator".
  // The conceptually easiest way is to create an MPI_Group that only
  // includes the shmem roots and then create a communicator from this
  // via MPI_Comm_create or MPI_Comm_create_group. The problem
  // with this is that we would have to exchange among all processes
  // which ones are shmem roots and which are not. This is awkward.
  //
  // A simpler way is to use MPI_Comm_split that uses "colors" to
  // indicate which sub-communicator each process wants to be in.
  // We use color=0 to indicate the group of shmem roots, and color=1
  // for all other processes -- the latter will simply not ever do
  // anything among themselves with the communicator so created.
  //
  // Using MPI_Comm_split has the additional benefit that, just as above,
  // we can choose where each rank will end up in shmem_roots_communicator.
  // We again set key=0 for the original root_process, and key=1 for all other
  // ranks; then, the global root becomes rank=0 on the
  // shmem_roots_communicator. We don't care how the other processes are
  // ordered.
  MPI_Comm shmem_roots_communicator;
  {
    const int key =
      (Utilities::MPI::this_mpi_process(communicator) == root_process ? 0 : 1);

    const int ierr = MPI_Comm_split(communicator,
                                    /*color=*/
                                    (is_shmem_root ? 0 : 1),
                                    key,
                                    &shmem_roots_communicator);
    AssertThrowMPI(ierr);

    // Again verify the explanation from above
    if (Utilities::MPI::this_mpi_process(communicator) == root_process)
      Assert(Utilities::MPI::this_mpi_process(shmem_roots_communicator) == 0,
             ExcInternalError());
  }

  const unsigned int shmem_roots_root_rank = 0;
  const bool         is_shmem_roots_root =
    (is_shmem_root && (Utilities::MPI::this_mpi_process(
                         shmem_roots_communicator) == shmem_roots_root_rank));

  // Now let the original root_process broadcast the current object to all
  // shmem roots. We know that the last rank is the original root process that
  // has all of the data.
  if (is_shmem_root)
    {
      if (std::is_trivial<T>::value)
        {
          // The data is "trivial", i.e., we can copy things directly without
          // having to go through the serialization/deserialization machinery of
          // Utilities::MPI::broadcast.
          //
          // In that case, first tell all of the other shmem roots how many
          // elements we will have to deal with, and let them resize their
          // (non-shared) arrays.
          const size_type new_size =
            Utilities::MPI::broadcast(shmem_roots_communicator,
                                      size(),
                                      shmem_roots_root_rank);
          if (is_shmem_roots_root == false)
            resize(new_size);

          // Then directly copy from the root process into these buffers
          int ierr = MPI_Bcast(elements.get(),
                               sizeof(T) * new_size,
                               MPI_CHAR,
                               shmem_roots_root_rank,
                               shmem_roots_communicator);
          AssertThrowMPI(ierr);
        }
      else
        {
          // The objects to be sent around are not "trivial", and so we have
          // to go through the serialization/deserialization machinery. On all
          // but the sending process, overwrite the current state with the
          // vector just broadcast.
          //
          // On the root rank, this would lead to resetting the 'entries'
          // pointer, which would trigger the deleter which would lead to a
          // deadlock. So we just send the result of the broadcast() call to
          // nirvana on the root process and keep our current state.
          if (Utilities::MPI::this_mpi_process(shmem_roots_communicator) == 0)
            Utilities::MPI::broadcast(shmem_roots_communicator,
                                      *this,
                                      shmem_roots_root_rank);
          else
            *this = Utilities::MPI::broadcast(shmem_roots_communicator,
                                              *this,
                                              shmem_roots_root_rank);
        }
    }

  // We no longer need the shmem roots communicator, so get rid of it
  {
    const int ierr = MPI_Comm_free(&shmem_roots_communicator);
    AssertThrowMPI(ierr);
  }


  // **** Step 3 ****
  // At this point, all shmem groups have one shmem root process that has
  // a copy of the data. This is the point where each shmem group should
  // establish a shmem area to put the data into. As mentioned above,
  // we know that the shmem roots are the last rank in their respective
  // shmem_group_communicator.
  //
  // The process for all of this works as follows: While all processes in
  // the shmem group participate in the generation of the shmem memory window,
  // only the shmem root actually allocates any memory -- the rest just
  // allocate zero bytes of their own. We allocate space for exactly
  // size() elements (computed on the shmem_root that already has the data)
  // and add however many bytes are necessary so that we know that we can align
  // things to 64-byte boundaries. The worst case happens if the memory system
  // gives us a pointer to an address one byte past a desired alignment
  // boundary, and in that case aligning the memory will require us to waste the
  // first (align_by-1) bytes. So we have to ask for
  //   size() * sizeof(T) + (align_by - 1)
  // bytes.
  //
  // Before MPI 4.0, there was no way to specify that we want memory aligned to
  // a certain number of bytes. This is going to come back to bite us further
  // down below when we try to get a properly aligned pointer to our memory
  // region, see the commentary there. Starting with MPI 4.0, one can set a
  // flag in an MPI_Info structure that requests a desired alignment, so we do
  // this for forward compatibility; MPI implementations ignore flags they don't
  // know anything about, and so setting this flag is backward compatible also
  // to older MPI versions.
  //
  // There is one final piece we can already take care of here. At the beginning
  // of all of this, only the shmem_root knows how many elements there are in
  // the array. But at the end of it, all processes of course need to know. We
  // could put this information somewhere into the shmem area, along with the
  // other data, but that seems clumsy. It turns out that when calling
  // MPI_Win_allocate_shared, we are asked for the value of a parameter called
  // 'disp_unit' whose meaning is difficult to determine from the MPI
  // documentation, and that we do not actually need. So we "abuse" it a bit: On
  // the shmem root, we put the array size into it. Later on, the remaining
  // processes can query the shmem root's value of 'disp_unit', and so will be
  // able to learn about the array size that way.
  MPI_Win        shmem_window;
  void *         base_ptr;
  const MPI_Aint align_by = 64;
  const MPI_Aint alloc_size =
    Utilities::MPI::broadcast(shmem_group_communicator,
                              (size() * sizeof(T) + (align_by - 1)),
                              0);

  {
    const int disp_unit = (is_shmem_root ? size() : 1);

    int ierr;

    MPI_Info mpi_info;
    ierr = MPI_Info_create(&mpi_info);
    AssertThrowMPI(ierr);
    ierr = MPI_Info_set(mpi_info,
                        "mpi_minimum_memory_alignment",
                        std::to_string(align_by).c_str());
    AssertThrowMPI(ierr);
    ierr = MPI_Win_allocate_shared((is_shmem_root ? alloc_size : 0),
                                   disp_unit,
                                   mpi_info,
                                   shmem_group_communicator,
                                   &base_ptr,
                                   &shmem_window);
    AssertThrowMPI(ierr);

    ierr = MPI_Info_free(&mpi_info);
    AssertThrowMPI(ierr);
  }


  // **** Step 4 ****
  // The next step is to teach all non-shmem root processes what the pointer to
  // the array is that the shmem-root created. MPI has a nifty way for this
  // given that only a single process actually allocated memory in the window:
  // When calling MPI_Win_shared_query, the MPI documentation says that
  // "When rank is MPI_PROC_NULL, the pointer, disp_unit, and size returned are
  // the pointer, disp_unit, and size of the memory segment belonging the lowest
  // rank that specified size > 0. If all processes in the group attached to the
  // window specified size = 0, then the call returns size = 0 and a baseptr as
  // if MPI_ALLOC_MEM was called with size = 0."
  //
  // This will allow us to obtain the pointer to the shmem root's memory area,
  // which is the only one we care about. (None of the other processes have
  // even allocated any memory.) But this will also retrieve the shmem root's
  // disp_unit, which in step 3 above we have abused to pass along the number of
  // elements in the array.
  //
  // We don't need to do this on the shmem root process: This process has
  // already gotten its base_ptr correctly set above, and we can determine the
  // array size by just calling size().
  unsigned int array_size =
    (is_shmem_root ? size() : numbers::invalid_unsigned_int);
  if (is_shmem_root == false)
    {
      int       disp_unit;
      MPI_Aint  alloc_size; // not actually used
      const int ierr = MPI_Win_shared_query(
        shmem_window, MPI_PROC_NULL, &alloc_size, &disp_unit, &base_ptr);
      AssertThrowMPI(ierr);

      // Make sure we actually got a pointer, and also unpack the array size as
      // discussed above.
      Assert(base_ptr != nullptr, ExcInternalError());

      array_size = disp_unit;
    }


  // **** Step 5 ****
  // Now that all processes know the address of the space that is visible to
  // everyone, we need to figure out whether it is properly aligned and if not,
  // find the next aligned address.
  //
  // std::align does that, but it also modifies its last two arguments. The
  // documentation of that function at
  // https://en.cppreference.com/w/cpp/memory/align is not entirely clear, but I
  // *think* that the following should do given that we do not use base_ptr and
  // available_space any further after the call to std::align.
  std::size_t available_space       = alloc_size;
  void *      base_ptr_backup       = base_ptr;
  T *         aligned_shmem_pointer = static_cast<T *>(
    std::align(align_by, array_size * sizeof(T), base_ptr, available_space));
  Assert(aligned_shmem_pointer != nullptr, ExcInternalError());

  // There is one step to guard against. It is *conceivable* that the base_ptr
  // we have previously obtained from MPI_Win_shared_query is mapped so
  // awkwardly into the different MPI processes' memory spaces that it is
  // aligned in one memory space, but not another. In that case, different
  // processes would align base_ptr differently, and adjust available_space
  // differently. We can check that by making sure that the max (or min) over
  // all processes is equal to every process's value. If that's not the case,
  // then the whole idea of aligning above is wrong and we need to rethink what
  // it means to align data in a shared memory space.
  //
  // One might be tempted to think that this is not how MPI implementations
  // actually arrange things. Alas, when developing this functionality in 2021,
  // this is really how at least OpenMPI ends up doing things. (This is with an
  // OpenMPI implementation of MPI 3.1, so it does not support the flag we set
  // in the MPI_Info structure above when allocating the memory window.) Indeed,
  // when running this code on three processes, one ends up with base_ptr values
  // of
  //     base_ptr=0x7f0842f02108
  //     base_ptr=0x7fc0a47881d0
  //     base_ptr=0x7f64872db108
  // which, most annoyingly, are aligned to 8 and 16 byte boundaries -- so there
  // is no common offset std::align could find that leads to a 64-byte
  // aligned memory address in all three memory spaces. That's a tremendous
  // nuisance and there is really nothing we can do about this other than just
  // fall back on the (unaligned) base_ptr in that case.
  if (Utilities::MPI::min(available_space, shmem_group_communicator) !=
      Utilities::MPI::max(available_space, shmem_group_communicator))
    aligned_shmem_pointer = static_cast<T *>(base_ptr_backup);


  // **** Step 6 ****
  // If this is the shmem root process, we need to copy the data into the
  // shared memory space.
  if (is_shmem_root)
    {
      if (std::is_trivial<T>::value == true)
        std::memcpy(aligned_shmem_pointer, elements.get(), sizeof(T) * size());
      else
        for (std::size_t i = 0; i < size(); ++i)
          new (&aligned_shmem_pointer[i]) T(std::move(elements[i]));
    }

  // Make sure that the shared memory host has copied the data before we try to
  // access it.
  MPI_Barrier(shmem_group_communicator);

  // **** Step 7 ****
  // Finally, we need to set the pointers of this object to what we just
  // learned. This also releases all memory that may have been in use
  // previously.
  //
  // The part that is a bit tricky is how to write the deleter of this
  // shared memory object. When we want to get rid of it, we need to
  // also release the MPI_Win object along with the shmem_group_communicator
  // object. That's because as long as we use the shared memory, we still need
  // to hold on to the MPI_Win object, and the MPI_Win object is based on the
  // communicator. (The former is definitely true, the latter is not quite clear
  // from the MPI documentation, but seems reasonable.) So we need to have a
  // deleter for the pointer that ensures that upon release of the memory, we
  // not only call the destructor of these memory elements (but only once, on
  // the shmem root!) but also destroy the MPI_Win and the communicator. All of
  // that is encapsulated in the following call where the deleter makes copies
  // of the arguments in the lambda capture.
  elements =
    decltype(elements)(aligned_shmem_pointer,
                       [is_shmem_root,
                        array_size,
                        aligned_shmem_pointer,
                        shmem_group_communicator,
                        shmem_window](T *) mutable {
                         if (is_shmem_root)
                           for (unsigned int i = 0; i < array_size; ++i)
                             aligned_shmem_pointer[i].~T();

                         int ierr;
                         ierr = MPI_Win_free(&shmem_window);
                         AssertThrowMPI(ierr);

                         ierr = MPI_Comm_free(&shmem_group_communicator);
                         AssertThrowMPI(ierr);
                       });

  // We then also have to set the other two pointers that define the state of
  // the current object. Note that the new buffer size is exactly as large as
  // necessary, i.e., can store size() elements, regardless of the number of
  // allocated elements in the original objects.
  used_elements_end      = elements.get() + array_size;
  allocated_elements_end = used_elements_end;

  // **** Consistency check ****
  // At this point, each process should have a copy of the data.
  // Verify this in some sort of round-about way
#      ifdef DEBUG
  const std::vector<char> packed_data = Utilities::pack(*this);
  const int               hash =
    std::accumulate(packed_data.begin(), packed_data.end(), int(0));
  Assert(Utilities::MPI::max(hash, communicator) == hash, ExcInternalError());
#      endif



#    else
  // If we only have MPI 2.x, then simply broadcast the current object to all
  // other processes and forego the idea of using shmem
  *this = Utilities::MPI::broadcast(communicator, *this, root_process);
#    endif
#  else
  // No MPI -> nothing to replicate
  (void)communicator;
  (void)root_process;
#  endif
}



template <class T>
inline void
AlignedVector<T>::swap(AlignedVector<T> &vec)
{
  std::swap(elements, vec.elements);
  std::swap(used_elements_end, vec.used_elements_end);
  std::swap(allocated_elements_end, vec.allocated_elements_end);
}



template <class T>
inline bool
AlignedVector<T>::empty() const
{
  return used_elements_end == elements.get();
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::size() const
{
  return used_elements_end - elements.get();
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::capacity() const
{
  return allocated_elements_end - elements.get();
}



template <class T>
inline typename AlignedVector<T>::reference AlignedVector<T>::
                                            operator[](const size_type index)
{
  AssertIndexRange(index, size());
  return elements[index];
}



template <class T>
inline typename AlignedVector<T>::const_reference AlignedVector<T>::
                                                  operator[](const size_type index) const
{
  AssertIndexRange(index, size());
  return elements[index];
}



template <typename T>
inline typename AlignedVector<T>::pointer
AlignedVector<T>::data()
{
  return elements.get();
}



template <typename T>
inline typename AlignedVector<T>::const_pointer
AlignedVector<T>::data() const
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::iterator
AlignedVector<T>::begin()
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::iterator
AlignedVector<T>::end()
{
  return used_elements_end;
}



template <class T>
inline typename AlignedVector<T>::const_iterator
AlignedVector<T>::begin() const
{
  return elements.get();
}



template <class T>
inline typename AlignedVector<T>::const_iterator
AlignedVector<T>::end() const
{
  return used_elements_end;
}



template <class T>
template <class Archive>
inline void
AlignedVector<T>::save(Archive &ar, const unsigned int) const
{
  size_type vec_size = size();
  ar &      vec_size;
  if (vec_size > 0)
    ar &boost::serialization::make_array(elements.get(), vec_size);
}



template <class T>
template <class Archive>
inline void
AlignedVector<T>::load(Archive &ar, const unsigned int)
{
  size_type vec_size = 0;
  ar &      vec_size;

  if (vec_size > 0)
    {
      reserve(vec_size);
      ar &boost::serialization::make_array(elements.get(), vec_size);
      used_elements_end = elements.get() + vec_size;
    }
}



template <class T>
inline typename AlignedVector<T>::size_type
AlignedVector<T>::memory_consumption() const
{
  size_type memory = sizeof(*this);
  for (const T *t = elements.get(); t != used_elements_end; ++t)
    memory += dealii::MemoryConsumption::memory_consumption(*t);
  memory += sizeof(T) * (allocated_elements_end - used_elements_end);
  return memory;
}


#endif // ifndef DOXYGEN


/**
 * Relational operator == for AlignedVector
 *
 * @relatesalso AlignedVector
 */
template <class T>
bool
operator==(const AlignedVector<T> &lhs, const AlignedVector<T> &rhs)
{
  if (lhs.size() != rhs.size())
    return false;
  for (typename AlignedVector<T>::const_iterator lit = lhs.begin(),
                                                 rit = rhs.begin();
       lit != lhs.end();
       ++lit, ++rit)
    if (*lit != *rit)
      return false;
  return true;
}



/**
 * Relational operator != for AlignedVector
 *
 * @relatesalso AlignedVector
 */
template <class T>
bool
operator!=(const AlignedVector<T> &lhs, const AlignedVector<T> &rhs)
{
  return !(operator==(lhs, rhs));
}


DEAL_II_NAMESPACE_CLOSE

#endif
