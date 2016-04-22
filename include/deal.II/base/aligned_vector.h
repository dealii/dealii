// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
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


#ifndef dealii__aligned_vector_h
#define dealii__aligned_vector_h

#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/type_traits.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parallel.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

#include <cstring>



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
template < class T >
class AlignedVector
{
public:
  /**
   * Declare standard types used in all containers. These types parallel those
   * in the <tt>C++</tt> standard libraries <tt>vector<...></tt> class.
   */
  typedef T                   value_type;
  typedef value_type         *pointer;
  typedef const value_type   *const_pointer;
  typedef value_type         *iterator;
  typedef const value_type   *const_iterator;
  typedef value_type         &reference;
  typedef const value_type   &const_reference;
  typedef std::size_t         size_type;

  /**
   * Empty constructor. Sets the vector size to zero.
   */
  AlignedVector ();

  /**
   * Sets the vector size to the given size and initializes all elements with
   * T().
   *
   * @dealiiOperationIsMultithreaded
   */
  AlignedVector (const size_type size,
                 const T        &init = T());

  /**
   * Destructor.
   */
  ~AlignedVector ();

  /**
   * Copy constructor.
   *
   * @dealiiOperationIsMultithreaded
   */
  AlignedVector (const AlignedVector<T> &vec);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move constructor. Create a new aligned vector by stealing the contents of
   * @p vec.
   *
   * @note This constructor is only available if deal.II is built with C++11.
   */
  AlignedVector (AlignedVector<T> &&vec);
#endif

  /**
   * Assignment to the input vector @p vec.
   *
   * @dealiiOperationIsMultithreaded
   */
  AlignedVector &
  operator = (const AlignedVector<T> &vec);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move assignment operator.
   *
   * @note This operator is only available if deal.II is built with C++11.
   */
  AlignedVector &
  operator = (AlignedVector<T> &&vec);
#endif

  /**
   * Change the size of the vector. It keeps old elements previously available
   * but does not initialize the newly allocated memory, leaving it in an
   * undefined state.
   *
   * @note This method can only be invoked for classes @p T that define a
   * default constructor, @p T(). Otherwise, compilation will fail.
   */
  void resize_fast (const size_type size);

  /**
   * Change the size of the vector. It keeps old elements previously
   * available, and initializes each element with the specified data. If the
   * new vector size is shorter than the old one, the memory is not released
   * unless the new size is zero.
   *
   * @dealiiOperationIsMultithreaded
   */
  void resize (const size_type size_in,
               const T        &init = T());

  /**
   * Reserve memory space for @p size elements. If the argument @p size is set
   * to zero, all previously allocated memory is released.
   *
   * In order to avoid too frequent reallocation (which involves copy of the
   * data), this function doubles the amount of memory occupied when the given
   * size is larger than the previously allocated size.
   */
  void reserve (const size_type size_alloc);

  /**
   * Releases all previously allocated memory and leaves the vector in a state
   * equivalent to the state after the default constructor has been called.
   */
  void clear ();

  /**
   * Inserts an element at the end of the vector, increasing the vector size
   * by one. Note that the allocated size will double whenever the previous
   * space is not enough to hold the new element.
   */
  void push_back (const T in_data);

  /**
   * Returns the last element of the vector (read and write access).
   */
  reference back ();

  /**
   * Returns the last element of the vector (read-only access).
   */
  const_reference back () const;

  /**
   * Inserts several elements at the end of the vector given by a range of
   * elements.
   */
  template <typename ForwardIterator>
  void insert_back (ForwardIterator begin,
                    ForwardIterator end);

  /**
   * Fills the vector with size() copies of the given input.
   *
   * @note This method can only be invoked for classes that define the copy
   * assignment operator. Otherwise, compilation will fail.
   *
   * @dealiiOperationIsMultithreaded
   */
  void fill (const T &element);

  /**
   * Swaps the given vector with the calling vector.
   */
  void swap (AlignedVector<T> &vec);

  /**
   * Returns whether the vector is empty, i.e., its size is zero.
   */
  bool empty () const;

  /**
   * Returns the size of the vector.
   */
  size_type size () const;

  /**
   * Returns the capacity of the vector, i.e., the size this vector can hold
   * without reallocation. Note that capacity() >= size().
   */
  size_type capacity () const;

  /**
   * Read-write access to entry @p index in the vector.
   */
  reference
  operator [] (const size_type index);

  /**
   * Read-only access to entry @p index in the vector.
   */
  const_reference operator [] (const size_type index) const;

  /**
   * Returns a read and write pointer to the beginning of the data array.
   */
  iterator begin ();

  /**
   * Returns a read and write pointer to the end of the data array.
   */
  iterator end ();

  /**
   * Returns a read-only pointer to the beginning of the data array.
   */
  const_iterator begin () const;

  /**
   * Returns a read-only pointer to the end of the data array.
   */
  const_iterator end () const;

  /**
   * Returns the memory consumption of the allocated memory in this class. If
   * the underlying type @p T allocates memory by itself, this memory is not
   * counted.
   */
  size_type memory_consumption () const;

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void save (Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void load (Archive &ar, const unsigned int version);

  BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

  /**
   * Pointer to actual class data.
   */
  T *_data;

  /**
   * Pointer to the end of valid data fields.
   */
  T *_end_data;

  /**
   * Pointer to the end of the allocated memory.
   */
  T *_end_allocated;
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
   * Move and class that actually issues the copy commands in AlignedVector.
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
   * @relates AlignedVector
   */
  template <typename T>
  class AlignedVectorMove : private parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size = 160000/sizeof(T)+1;
  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise works in serial. Copies the data from the half-open
     * interval between @p source_begin and @p source_end to array starting at
     * @p destination (by calling the copy constructor with placement new). If
     * the flag copy_source is set to @p true, the elements from the source
     * array are simply copied. If it is set to @p false, the data is moved
     * between the two arrays by invoking the destructor on the source range
     * (preparing for a subsequent call to free).
     */
    AlignedVectorMove (T *source_begin,
                       T *source_end,
                       T *destination,
                       const bool copy_source)
      :
      source_ (source_begin),
      destination_ (destination),
      copy_source_ (copy_source)
    {
      Assert (source_end >= source_begin, ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        apply_to_subrange (0, size);
      else
        apply_parallel (0, size, minimum_parallel_grain_size);
    }

    /**
     * This method moves elements from the source to the destination given in
     * the constructor on a subrange given by two integers.
     */
    virtual void apply_to_subrange (const std::size_t begin,
                                    const std::size_t end) const
    {
      // for classes trivial assignment can use memcpy. cast element to
      // (void*) to silence compiler warning for virtual classes (they will
      // never arrive here because they are non-trivial).

      if (std_cxx11::is_trivial<T>::value == true)
        std::memcpy ((void *)(destination_+begin), source_+begin,
                     (end-begin)*sizeof(T));
      else if (copy_source_ == false)
        for (std::size_t i=begin; i<end; ++i)
          {
            // initialize memory (copy construct by placement new), and
            // destruct the source
            new (&destination_[i]) T(source_[i]);
            source_[i].~T();
          }
      else
        for (std::size_t i=begin; i<end; ++i)
          new (&destination_[i]) T(source_[i]);
    }

  private:
    T *source_;
    T *destination_;
    const bool copy_source_;
  };

  /**
   * Class that issues the set commands for AlignedVector.
   *
   * @tparam initialize_memory Sets whether the the set command should
   * initialize memory (with a call to the copy constructor) or rather use the
   * copy assignment operator. A template is necessary to select the
   * appropriate operation since some classes might define only one of those
   * two operations.
   *
   * @relates AlignedVector
   */
  template <typename T, bool initialize_memory>
  class AlignedVectorSet : private parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size = 160000/sizeof(T)+1;
  public:
    /**
     * Constructor. Issues a parallel call if there are sufficiently many
     * elements, otherwise work in serial.
     */
    AlignedVectorSet (const std::size_t size,
                      const T &element,
                      T *destination)
      :
      element_ (element),
      destination_ (destination),
      trivial_element (false)
    {
      if (size == 0)
        return;

      // do not use memcmp for long double because on some systems it does not
      // completely fill its memory and may lead to false positives in
      // e.g. valgrind
      if (std_cxx11::is_trivial<T>::value == true &&
          types_are_equal<T,long double>::value == false)
        {
          const unsigned char zero [sizeof(T)] = {};
          // cast element to (void*) to silence compiler warning for virtual
          // classes (they will never arrive here because they are
          // non-trivial).
          if (std::memcmp(zero, (void *)&element, sizeof(T)) == 0)
            trivial_element = true;
        }
      if (size < minimum_parallel_grain_size)
        apply_to_subrange (0, size);
      else
        apply_parallel (0, size, minimum_parallel_grain_size);
    }

    /**
     * This sets elements on a subrange given by two integers.
     */
    virtual void apply_to_subrange (const std::size_t begin,
                                    const std::size_t end) const
    {
      // for classes with trivial assignment of zero can use memset. cast
      // element to (void*) to silence compiler warning for virtual
      // classes (they will never arrive here because they are
      // non-trivial).
      if (std_cxx11::is_trivial<T>::value == true && trivial_element)
        std::memset ((void *)(destination_+begin), 0, (end-begin)*sizeof(T));
      else
        copy_construct_or_assign(begin, end,
                                 ::dealii::internal::bool2type<initialize_memory>());
    }

  private:
    const T &element_;
    mutable T *destination_;
    bool trivial_element;

    // copy assignment operation
    void copy_construct_or_assign(const std::size_t begin,
                                  const std::size_t end,
                                  ::dealii::internal::bool2type<false>) const
    {
      for (std::size_t i=begin; i<end; ++i)
        destination_[i] = element_;
    }

    // copy constructor (memory initialization)
    void copy_construct_or_assign(const std::size_t begin,
                                  const std::size_t end,
                                  ::dealii::internal::bool2type<true>) const
    {
      for (std::size_t i=begin; i<end; ++i)
        new (&destination_[i]) T(element_);
    }
  };

} // end of namespace internal


#ifndef DOXYGEN


template < class T >
inline
AlignedVector<T>::AlignedVector ()
  :
  _data (0),
  _end_data (0),
  _end_allocated (0)
{}



template < class T >
inline
AlignedVector<T>::AlignedVector (const size_type size,
                                 const T &init)
  :
  _data (0),
  _end_data (0),
  _end_allocated (0)
{
  if (size > 0)
    resize (size, init);
}



template < class T >
inline
AlignedVector<T>::~AlignedVector ()
{
  clear();
}



template < class T >
inline
AlignedVector<T>::AlignedVector (const AlignedVector<T> &vec)
  :
  _data (0),
  _end_data (0),
  _end_allocated (0)
{
  // copy the data from vec
  reserve (vec._end_data - vec._data);
  _end_data = _end_allocated;
  internal::AlignedVectorMove<T> (vec._data, vec._end_data, _data, true);
}



#ifdef DEAL_II_WITH_CXX11
template < class T >
inline
AlignedVector<T>::AlignedVector (AlignedVector<T> &&vec)
  :
  _data (vec._data),
  _end_data (vec._end_data),
  _end_allocated (vec._end_allocated)
{
  vec._data = nullptr;
  vec._end_data = nullptr;
  vec._end_allocated = nullptr;
}
#endif



template < class T >
inline
AlignedVector<T> &
AlignedVector<T>::operator = (const AlignedVector<T> &vec)
{
  resize(0);
  resize_fast (vec._end_data - vec._data);
  internal::AlignedVectorMove<T> (vec._data, vec._end_data, _data, true);
  return *this;
}



#ifdef DEAL_II_WITH_CXX11
template < class T >
inline
AlignedVector<T> &
AlignedVector<T>::operator = (AlignedVector<T> &&vec)
{
  clear();

  _data = vec._data;
  _end_data = vec._end_data;
  _end_allocated = vec._end_allocated;

  vec._data = nullptr;
  vec._end_data = nullptr;
  vec._end_allocated = nullptr;

  return *this;
}
#endif



template < class T >
inline
void
AlignedVector<T>::resize_fast (const size_type size_in)
{
  const size_type old_size = size();
  if (std_cxx11::is_trivial<T>::value == false && size_in < old_size)
    {
      // call destructor on fields that are released. doing it backward
      // releases the elements in reverse order as compared to how they were
      // created
      while (_end_data != _data+size_in)
        (--_end_data)->~T();
    }
  reserve (size_in);
  _end_data = _data + size_in;

  // need to still set the values in case the class is non-trivial because
  // virtual classes etc. need to run their (default) constructor
  if (std_cxx11::is_trivial<T>::value == false && size_in > old_size)
    dealii::internal::AlignedVectorSet<T,true> (size_in-old_size, T(), _data+old_size);
}



template < class T >
inline
void
AlignedVector<T>::resize (const size_type size_in,
                          const T        &init)
{
  const size_type old_size = size();
  if (std_cxx11::is_trivial<T>::value == false && size_in < old_size)
    {
      // call destructor on fields that are released. doing it backward
      // releases the elements in reverse order as compared to how they were
      // created
      while (_end_data != _data+size_in)
        (--_end_data)->~T();
    }
  reserve (size_in);
  _end_data = _data + size_in;

  // finally set the desired init values
  if (size_in > old_size)
    dealii::internal::AlignedVectorSet<T,true> (size_in-old_size, init, _data+old_size);
}



template < class T >
inline
void
AlignedVector<T>::reserve (const size_type size_alloc)
{
  const size_type old_size = _end_data - _data;
  const size_type allocated_size = _end_allocated - _data;
  if (size_alloc > allocated_size)
    {
      // if we continuously increase the size of the vector, we might be
      // reallocating a lot of times. therefore, try to increase the size more
      // aggressively
      size_type new_size = size_alloc;
      if (size_alloc < (2 * allocated_size))
        new_size = 2 * allocated_size;

      const size_type size_actual_allocate = new_size * sizeof(T);

      // allocate and align along 64-byte boundaries (this is enough for all
      // levels of vectorization currently supported by deal.II)
      T *new_data;
      Utilities::System::posix_memalign ((void **)&new_data, 64, size_actual_allocate);

      // copy data in case there was some content before and release the old
      // memory with the function corresponding to the one used for allocating
      std::swap (_data, new_data);
      _end_data = _data + old_size;
      _end_allocated = _data + new_size;
      if (_end_data != _data)
        {
          dealii::internal::AlignedVectorMove<T>(new_data, new_data + old_size,
                                                 _data, false);
          free(new_data);
        }
      else
        Assert(new_data == 0, ExcInternalError());
    }
  else if (size_alloc == 0)
    clear();
}



template < class T >
inline
void
AlignedVector<T>::clear ()
{
  if (_data != 0)
    {
      if (std_cxx11::is_trivial<T>::value == false)
        while (_end_data != _data)
          (--_end_data)->~T();

      free(_data);
    }
  _data = 0;
  _end_data = 0;
  _end_allocated = 0;
}



template < class T >
inline
void
AlignedVector<T>::push_back (const T in_data)
{
  Assert (_end_data <= _end_allocated, ExcInternalError());
  if (_end_data == _end_allocated)
    reserve (std::max(2*capacity(),static_cast<size_type>(16)));
  if (std_cxx11::is_trivial<T>::value == false)
    new (_end_data) T;
  *_end_data++ = in_data;
}



template < class T >
inline
typename AlignedVector<T>::reference
AlignedVector<T>::back ()
{
  AssertIndexRange (0, size());
  T *field = _end_data - 1;
  return *field;
}



template < class T >
inline
typename AlignedVector<T>::const_reference
AlignedVector<T>::back () const
{
  AssertIndexRange (0, size());
  const T *field = _end_data - 1;
  return *field;
}



template < class T >
template <typename ForwardIterator>
inline
void
AlignedVector<T>::insert_back (ForwardIterator begin,
                               ForwardIterator end)
{
  const unsigned int old_size = size();
  reserve (old_size + (end-begin));
  for ( ; begin != end; ++begin, ++_end_data)
    {
      if (std_cxx11::is_trivial<T>::value == false)
        new (_end_data) T;
      *_end_data = *begin;
    }
}



template < class T >
inline
void
AlignedVector<T>::fill (const T &value)
{
  dealii::internal::AlignedVectorSet<T,false> (size(), value, _data);
}



template < class T >
inline
void
AlignedVector<T>::swap (AlignedVector<T> &vec)
{
  std::swap (_data, vec._data);
  std::swap (_end_data, vec._end_data);
  std::swap (_end_allocated, vec._end_allocated);
}



template < class T >
inline
bool
AlignedVector<T>::empty () const
{
  return _end_data == _data;
}



template < class T >
inline
typename AlignedVector<T>::size_type
AlignedVector<T>::size () const
{
  return _end_data - _data;
}



template < class T >
inline
typename AlignedVector<T>::size_type
AlignedVector<T>::capacity () const
{
  return _end_allocated - _data;
}



template < class T >
inline
typename AlignedVector<T>::reference
AlignedVector<T>::operator [] (const size_type index)
{
  AssertIndexRange (index, size());
  return _data[index];
}



template < class T >
inline
typename AlignedVector<T>::const_reference
AlignedVector<T>::operator [] (const size_type index) const
{
  AssertIndexRange (index, size());
  return _data[index];
}



template < class T >
inline
typename AlignedVector<T>::iterator
AlignedVector<T>::begin ()
{
  return _data;
}



template < class T >
inline
typename AlignedVector<T>::iterator
AlignedVector<T>::end ()
{
  return _end_data;
}



template < class T >
inline
typename AlignedVector<T>::const_iterator
AlignedVector<T>::begin () const
{
  return _data;
}



template < class T >
inline
typename AlignedVector<T>::const_iterator
AlignedVector<T>::end () const
{
  return _end_data;
}



template < class T >
template < class Archive >
inline
void
AlignedVector<T>::save (Archive &ar, const unsigned int) const
{
  size_type vec_size (size());
  ar &vec_size;
  if (vec_size > 0)
    ar &boost::serialization::make_array(_data, vec_size);
}



template < class T >
template < class Archive >
inline
void
AlignedVector<T>::load (Archive &ar, const unsigned int)
{
  size_type vec_size = 0;
  ar &vec_size ;

  if (vec_size > 0)
    {
      reserve(vec_size);
      ar &boost::serialization::make_array(_data, vec_size);
      _end_data = _data + vec_size;
    }
}



template < class T >
inline
typename AlignedVector<T>::size_type
AlignedVector<T>::memory_consumption () const
{
  size_type memory = sizeof(*this);
  for (const T *t = _data ; t != _end_data; ++t)
    memory += dealii::MemoryConsumption::memory_consumption(*t);
  memory += sizeof(T) * (_end_allocated-_end_data);
  return memory;
}


#endif // ifndef DOXYGEN


/**
 * Relational operator == for AlignedVector
 *
 * @relates AlignedVector
 */
template < class T >
bool operator == (const AlignedVector<T> &lhs,
                  const AlignedVector<T> &rhs)
{
  if (lhs.size() != rhs.size())
    return false;
  for (typename AlignedVector<T>::const_iterator lit = lhs.begin(),
       rit = rhs.begin(); lit != lhs.end(); ++lit, ++rit)
    if (*lit != *rit)
      return false;
  return true;
}




/**
 * Relational operator != for AlignedVector
 *
 * @relates AlignedVector
 */
template < class T >
bool operator != (const AlignedVector<T> &lhs,
                  const AlignedVector<T> &rhs)
{
  return !(operator==(lhs, rhs));
}


DEAL_II_NAMESPACE_CLOSE

#endif
