// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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

#ifndef dealii__memory_consumption_h
#define dealii__memory_consumption_h


#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/std_cxx11/type_traits.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/base/std_cxx11/array.h>

#include <string>
#include <complex>
#include <vector>
#include <cstddef>
#include <cstring>

DEAL_II_NAMESPACE_OPEN


// forward declaration
template <typename T> class VectorizedArray;


/**
 * This namespace provides functions helping to determine the amount of memory
 * used by objects. The goal is not necessarily to give the amount of memory
 * used up to the last bit (what is the memory used by a <tt>std::map</tt>
 * object?), but rather to aid in the search for memory bottlenecks.
 *
 * This namespace has a single member function memory_consumption() and a lot
 * of specializations. Depending on the argument type of the function, there
 * are several modes of operation:
 *
 * <ol>
 * <li> If the argument is a fundamental C++ data type (such as <tt>bool</tt>,
 * <tt>float</tt>, <tt>double</tt> or any of the integer types), then
 * memory_consumption() just returns <tt>sizeof</tt> of its argument. The
 * library also provides an estimate for the amount of memory occupied by a
 * <tt>std::string</tt>.
 *
 * <li> For objects, which are neither standard types, nor vectors,
 * memory_consumption() will simply call the member function of same name. It
 * is up to the implementation of the data type to provide a good estimate of
 * the amount of memory used. Inside this function, the use of
 * MemoryConsumption::memory_consumption() for compounds of the class helps to
 * obtain this estimate. Most classes in the deal.II library have such a
 * member function.
 *
 * <li> For vectors and C++ arrays of objects, memory_consumption()
 * recursively calls itself for all entries and adds the results to the size
 * of the object itself. Some optimized specializations for standard data
 * types exist.
 *
 * <li> For vectors of regular pointers, memory_consumption(T*) returns the
 * size of the vector of pointers, ignoring the size of the objects.
 *
 * </ol>
 *
 * <h3>Extending this namespace</h3>
 *
 * The function in this namespace and the functionality provided by it relies
 * on the assumption that there is either a function
 * <tt>memory_consumption(T)</tt> in this namespace determining the amount of
 * memory used by objects of type <tt>T</tt> or that the class <tt>T</tt> has
 * a member function of that name. While the latter is true for almost all
 * classes in deal.II, we have only implemented the first kind of functions
 * for the most common data types, such as fundamental types, strings, C++
 * vectors, C-style arrays, and C++ pairs. These functions therefore do not
 * cover, for example, C++ maps, lists, etc. If you need such functions feel
 * free to implement them and send them to us for inclusion.
 *
 * @ingroup memory
 * @author Wolfgang Bangerth, documentation updated by Guido Kanschat, David
 * Wells
 * @date 2000, 2015
 */
namespace MemoryConsumption
{
  /**
   * Calculate the memory consumption of a fundamental type. See
   * EnableIfScalar for a discussion on how this restriction (SFINAE) is
   * implemented.
   */
  template <typename T>
  inline
  typename std_cxx11::enable_if<std_cxx11::is_fundamental<T>::value, std::size_t>::type
  memory_consumption (const T &t);

  /**
   * Estimate the memory consumption of an object. If no further template
   * specialization (past this one) is available for the type <tt>T</tt>, then
   * this function returns the member function
   * <tt>t.memory_consumption()</tt>'s value.
   */
  template <typename T>
  inline
  typename std_cxx11::enable_if<!(std_cxx11::is_fundamental<T>::value || std_cxx11::is_pointer<T>::value), std::size_t>::type
  memory_consumption (const T &t);

  /**
   * Determine the amount of memory consumed by a C-style string. The returned
   * value does not include the size of the pointer. This function only
   * measures up to (and including) the NUL byte; the underlying buffer may be
   * larger.
   */
  inline
  std::size_t memory_consumption (const char *string);

  /**
   * Determine the amount of memory in bytes consumed by a
   * <tt>std::complex</tt> variable.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std::complex<T> &);

  /**
   * Determine the amount of memory in bytes consumed by a
   * <tt>VectorizedArray</tt> variable.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const VectorizedArray<T> &);

  /**
   * Determine an estimate of the amount of memory in bytes consumed by a
   * <tt>std::string</tt> variable.
   */
  inline
  std::size_t memory_consumption (const std::string &s);

  /**
   * Determine the amount of memory in bytes consumed by a
   * <tt>std::vector</tt> of elements of type <tt>T</tt> by
   * calling memory_consumption() for each entry.
   *
   * This function loops over all entries of the vector and determines their
   * sizes using memory_consumption() for each <tt>v[i]</tt>. If the entries
   * are of constant size, there might be another global function
   * memory_consumption() for this data type or if there is a member function
   * of that class of that names that returns a constant value and the
   * compiler will unroll this loop so that the operation is fast. If the size
   * of the data elements is variable, for example if they do memory
   * allocation themselves, then the operation will necessarily be more
   * expensive.
   *
   * Using the algorithm, in particular the loop over all elements, it is
   * possible to also compute the memory consumption of vectors of vectors,
   * vectors of strings, etc, where the individual elements may have vastly
   * different sizes.
   *
   * Note that this algorithm also takes into account the size of elements
   * that are allocated by this vector but not currently used.
   *
   * For the most commonly used vectors, there are special functions that
   * compute the size without a loop. This also applies for the special case
   * of vectors of bools.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std::vector<T> &v);

  /**
   * Determine the amount of memory in bytes consumed by a
   * <tt>std_cxx11::array</tt> of <tt>N</tt> elements of type <tt>T</tt> by
   * calling memory_consumption() for each entry.
   *
   * This function loops over all entries of the array and determines their
   * sizes using memory_consumption() for each <tt>v[i]</tt>. If the entries
   * are of constant size, there might be another global function
   * memory_consumption() for this data type or if there is a member function
   * of that class of that names that returns a constant value and the
   * compiler will unroll this loop so that the operation is fast. If the size
   * of the data elements is variable, for example if they do memory
   * allocation themselves, then the operation will necessarily be more
   * expensive.
   *
   * Using the algorithm, in particular the loop over all elements, it is
   * possible to also compute the memory consumption of arrays of vectors,
   * arrays of strings, etc, where the individual elements may have vastly
   * different sizes.
   */
  template <typename T, int N>
  inline
  std::size_t memory_consumption (const std_cxx11::array<T,N> &v);

  /**
   * Estimate the amount of memory (in bytes) occupied by a C-style array.
   * Since in this library we do not usually store simple data elements like
   * <tt>double</tt>s in such arrays (but rather use <tt>std::vector</tt>s or
   * deal.II <tt>Vector</tt> objects), we do not provide specializations like
   * for the <tt>std::vector</tt> arrays, but always use the loop over all
   * elements.
   */
  template <typename T, int N>
  inline
  std::size_t memory_consumption (const T (&v)[N]);

  /**
   * Specialization of the determination of the memory consumption of a
   * vector, here for a vector of <tt>bool</tt>s.
   *
   * This is a special case, as the bools are not stored one-by-one, but as a
   * bit field.
   */
  inline
  std::size_t memory_consumption (const std::vector<bool> &v);

  /**
   * Determine an estimate of the amount of memory in bytes consumed by a pair
   * of values.
   */
  template <typename A, typename B>
  inline
  std::size_t memory_consumption (const std::pair<A,B> &p);

  /**
   * Calculate the memory consumption of a pointer.
   *
   * @note This function is overloaded for C-style strings; see the
   * documentation of that function for that case.
   *
   * @note This returns the size of the pointer, not the size of the object
   * pointed to.
   */
  template<typename T>
  inline
  std::size_t
  memory_consumption (const T *const);

  /**
   * Return the amount of memory used by a shared pointer.
   *
   * @note This returns the size of the pointer, not of the object pointed to.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std_cxx11::shared_ptr<T> &);

  /**
   * Return the amount of memory used by a std_cxx11::unique_ptr object.
   *
   * @note This returns the size of the pointer, not of the object pointed to.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std_cxx11::unique_ptr<T> &);
}



// now comes the implementation of these functions

namespace MemoryConsumption
{
  template <typename T>
  inline
  typename std_cxx11::enable_if<std_cxx11::is_fundamental<T>::value, std::size_t>::type
  memory_consumption(const T &)
  {
    return sizeof(T);
  }



  inline
  std::size_t memory_consumption (const char *string)
  {
    if (string == NULL)
      {
        return 0;
      }
    else
      {
        return sizeof(char)*(strlen(string) /*Remember the NUL*/ + 1);
      }
  }



  template <typename T>
  inline
  std::size_t memory_consumption (const std::complex<T> &)
  {
    return sizeof(std::complex<T>);
  }



  template <typename T>
  inline
  std::size_t memory_consumption (const VectorizedArray<T> &)
  {
    return sizeof(VectorizedArray<T>);
  }



  inline
  std::size_t memory_consumption (const std::string &s)
  {
    return sizeof(s) + s.length();
  }



  template <typename T>
  std::size_t memory_consumption (const std::vector<T> &v)
  {
    // shortcut for types that do not allocate memory themselves
    if (std_cxx11::is_fundamental<T>::value || std_cxx11::is_pointer<T>::value)
      {
        return v.capacity()*sizeof(T) + sizeof(v);
      }
    else
      {
        std::size_t mem = sizeof(std::vector<T>);
        for (unsigned int i=0; i<v.size(); ++i)
          {
            mem += memory_consumption(v[i]);
          }
        mem += (v.capacity() - v.size())*sizeof(T);
        return mem;
      }
  }



  template <typename T, std::size_t N>
  std::size_t memory_consumption (const std_cxx11::array<T,N> &v)
  {
    // shortcut for types that do not allocate memory themselves
    if (std_cxx11::is_fundamental<T>::value || std_cxx11::is_pointer<T>::value)
      {
        return sizeof(v);
      }
    else
      {
        std::size_t mem = 0;
        for (std::size_t i=0; i!=N; ++i)
          mem += memory_consumption(v[i]);
        return mem;
      }
  }



  template <typename T, int N>
  std::size_t memory_consumption (const T (&v)[N])
  {
    std::size_t mem = 0;
    for (unsigned int i=0; i<N; ++i)
      mem += memory_consumption(v[i]);
    return mem;
  }



  inline
  std::size_t memory_consumption (const std::vector<bool> &v)
  {
    return v.capacity() / 8 + sizeof(v);
  }



  template <typename A, typename B>
  inline
  std::size_t memory_consumption (const std::pair<A,B> &p)
  {
    return (memory_consumption(p.first) +
            memory_consumption(p.second));
  }



  template<typename T>
  inline
  std::size_t
  memory_consumption(const T *const)
  {
    return sizeof(T *);
  }



  template <typename T>
  inline
  std::size_t
  memory_consumption (const std_cxx11::shared_ptr<T> &)
  {
    return sizeof(std_cxx11::shared_ptr<T>);
  }



  template <typename T>
  inline
  std::size_t
  memory_consumption (const std_cxx11::unique_ptr<T> &)
  {
    return sizeof(std_cxx11::unique_ptr<T>);
  }



  template <typename T>
  inline
  typename std_cxx11::enable_if<!(std_cxx11::is_fundamental<T>::value || std_cxx11::is_pointer<T>::value), std::size_t>::type
  memory_consumption (const T &t)
  {
    return t.memory_consumption();
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
