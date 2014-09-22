// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#ifndef __deal2__memory_consumption_h
#define __deal2__memory_consumption_h


#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <string>
#include <complex>
#include <vector>
#include <cstddef>

DEAL_II_NAMESPACE_OPEN


// forward declaration
template <typename T> class VectorizedArray;


/**
 * This namespace provides functions helping to determine the amount
 * of memory used by objects. The goal is not necessarily to give the
 * amount of memory used up to the last bit (what is the memory used
 * by an STL <tt>std::map<></tt> object?), but rather to aid in the search for
 * memory bottlenecks.
 *
 * This namespace has a single member function memory_consumption()
 * and a lot of specializations. Depending on the argument type of the
 * function, there are several modes of operation:
 *
 * <ol>
 * <li> The argument is a standard C++ data type, namely,
 * <tt>bool</tt>, <tt>float</tt>, <tt>double</tt> or any of the
 * integer types. In that case, memory_consumption() simple returns
 * <tt>sizeof</tt> of its argument. The library also provides an
 * estimate for the amount of memory occupied by a
 * <tt>std::string</tt> this way.
 *
 * <li> For objects, which are neither standard types, nor vectors,
 * memory_consumption() will simply call the member function of same
 * name. It is up to the implementation of the data type to provide a
 * good estimate of the amount of memory used. Inside this function,
 * the use of MemoryConsumpton::memory_consumption() for compounds of
 * the class helps to obtain this estimate. Most classes in the
 * deal.II library have such a member function.
 *
 * <li> For vectors and C++ arrays of objects, memory_consumption()
 * recursively calls itself for all entries and adds the results to
 * the size of the object itself. Some optimized specializations for
 * standard data types exist.
 *
 * <li> For vectors of regular pointers, memory_consumption(T*)
 * returns the size of the vector of pointers, ignoring the size of
 * the objects.
 *
 * </ol>
 *
 * <h3>Extending this namespace</h3>
 *
 * The function in this namespace and the functionality provided by
 * it relies on the assumption that there is either a specialized function
 * <tt>memory_consumption(T)</tt> in this namespace determining the amount
 * of memory used by objects of type <tt>T</tt>, or that the class <tt>T</tt> has
 * a  member function of that name. While the latter is
 * true for almost all classes in deal.II, we have only implemented the
 * first kind of functions for the most common data types, such as
 * atomic types, strings, C++ vectors, C-style arrays, and C++
 * pairs. These functions therefore do not cover, for example, C++
 * maps, lists, etc. If you need such functions feel free to implement
 * them and send them to us for inclusion.
 *
 * @ingroup memory
 * @author Wolfgang Bangerth, documentation updated by Guido Kanschat
 * @date 2000
 */
namespace MemoryConsumption
{

  /**
   * This function is the generic
   * interface for determining the
   * memory used by an object. If no
   * specialization for the type
   * <tt>T</tt> is specified, it will
   * call the member function
   * <tt>t.memory_consumption()</tt>.
   *
   * The library provides
   * specializations for all basic
   * C++ data types. Every additional
   * type needs to have a member
   * function memory_consumption()
   * callable for constant objects to
   * be used in this framework.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const T &t);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>bool</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const bool);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>char</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const char);

  /**
   * Determine the amount of memory
   * in bytes consumed by a
   * <tt>short int</tt> variable.
   */
  inline
  std::size_t memory_consumption (const short int);

  /**
   * Determine the amount of memory
   * in bytes consumed by a
   * <tt>short unsigned int</tt> variable.
   */
  inline
  std::size_t memory_consumption (const short unsigned int);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>int</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const int);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>unsigned int</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const unsigned int);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>unsigned long long int</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const unsigned long long int);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>float</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const float);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>double</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const double);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>long double</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const long double);

  /**
   * Determine the amount of memory
   * in bytes consumed by a <tt>std::complex</tt>
   * variable.
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
   * Determine an estimate of the
   * amount of memory in bytes
   * consumed by a <tt>std::string</tt>
   * variable.
   */
  inline
  std::size_t memory_consumption (const std::string &s);

  /**
   * Determine the amount of memory
   * in bytes consumed by a
   * <tt>std::vector</tt> of elements
   * of type <tt>T</tt> by
   * recursively calling
   * memory_consumption() for each entry.
   *
   * This function loops over all
   * entries of the vector and
   * determines their sizes using
   * memory_consumption() for each
   * <tt>v[i]</tt>. If the entries
   * are of constant size, there
   * might be another global function
   * memory_consumption() for this
   * data type or if there is a
   * member function of that class of
   * that names that returns a
   * constant value and the compiler
   * will unroll this loop so that
   * the operation is fast. If the
   * size of the data elements is
   * variable, for example if they do
   * memory allocation themselves,
   * then the operation will
   * necessarily be more expensive.
   *
   * Using the algorithm, in
   * particular the loop over all
   * elements, it is possible to also
   * compute the memory consumption
   * of vectors of vectors, vectors
   * of strings, etc, where the
   * individual elements may have
   * vastly different sizes.
   *
   * Note that this algorithm also
   * takes into account the size of
   * elements that are allocated by
   * this vector but not currently
   * used.
   *
   * For the most commonly used
   * vectors, there are special
   * functions that compute the size
   * without a loop. This also
   * applies for the special case of
   * vectors of bools.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std::vector<T> &v);

  /**
   * Estimate the amount of memory
   * (in bytes) occupied by a
   * C-style array. Since in this
   * library we do not usually
   * store simple data elements
   * like <tt>double</tt>s in such
   * arrays (but rather use STL
   * <tt>std::vector</tt>s or deal.II
   * <tt>Vector</tt> objects), we do not
   * provide specializations like
   * for the <tt>std::vector</tt> arrays, but
   * always use the loop over all
   * elements.
   */
  template <typename T, int N>
  inline
  std::size_t memory_consumption (const T (&v)[N]);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>bool</tt>s.
   *
   * This is a special case, as the
   * bools are not stored
   * one-by-one, but as a bit
   * field.
   */
  inline
  std::size_t memory_consumption (const std::vector<bool> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>int</tt>s.
   */
  inline
  std::size_t memory_consumption (const std::vector<int> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>double</tt>s.
   */
  inline
  std::size_t memory_consumption (const std::vector<double> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>float</tt>s.
   */
  inline
  std::size_t memory_consumption (const std::vector<float> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>char</tt>s.
   */
  inline
  std::size_t memory_consumption (const std::vector<char> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of <tt>unsigned char</tt>s.
   */
  inline
  std::size_t memory_consumption (const std::vector<unsigned char> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of pointers.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std::vector<T *> &v);

  /**
   * Specialization of the
   * determination of the memory
   * consumption of a vector, here
   * for a vector of strings. This
   * function is not necessary from a
   * strict C++ viewpoint, since it
   * could be generated, but is
   * necessary for compatibility with
   * IBM's xlC 5.0 compiler, and
   * doesn't harm for other compilers
   * as well.
   */
  std::size_t memory_consumption (const std::vector<std::string> &v);


  /**
   * Determine an estimate of the
   * amount of memory in bytes
   * consumed by a pair of values.
   */
  template <typename A, typename B>
  inline
  std::size_t memory_consumption (const std::pair<A,B> &p);

  /**
   * Return the amount of memory
   * used by a pointer.
   *
   * @note This returns the size of
   * the pointer, not of the object
   * pointed to.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const T *const);

  /**
   * Return the amount of memory
   * used by a pointer.
   *
   * @note This returns the size of
   * the pointer, not of the object
   * pointed to.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (T *const);

  /**
   * Return the amount of memory
   * used by a void pointer.
   *
   * Note that we needed this
   * function since <tt>void</tt> is no
   * type and a <tt>void*</tt> is thus
   * not caught by the general
   * <tt>T*</tt> template function
   * above.
   *
   * @note This returns the size of
   * the pointer, not of the object
   * pointed to.
   */
  inline
  std::size_t memory_consumption (void *const);

  /**
   * Return the amount of memory used
   * by a shared pointer.
   *
   * @note This returns the size of
   * the pointer, not of the object
   * pointed to.
   */
  template <typename T>
  inline
  std::size_t memory_consumption (const std_cxx11::shared_ptr<T> &);
}



// now comes the implementation of these functions

namespace MemoryConsumption
{
  inline
  std::size_t memory_consumption (const bool)
  {
    return sizeof(bool);
  }



  inline
  std::size_t memory_consumption (const char)
  {
    return sizeof(char);
  }



  inline
  std::size_t memory_consumption (const short int)
  {
    return sizeof(short int);
  }



  inline
  std::size_t memory_consumption (const short unsigned int)
  {
    return sizeof(short unsigned int);
  }



  inline
  std::size_t memory_consumption (const int)
  {
    return sizeof(int);
  }



  inline
  std::size_t memory_consumption (const unsigned int)
  {
    return sizeof(unsigned int);
  }



  inline
  std::size_t memory_consumption (const unsigned long int)
  {
    return sizeof(unsigned long int);
  }



  inline
  std::size_t memory_consumption (const unsigned long long int)
  {
    return sizeof(unsigned long long int);
  }



  inline
  std::size_t memory_consumption (const float)
  {
    return sizeof(float);
  }



  inline
  std::size_t memory_consumption (const double)
  {
    return sizeof(double);
  }



  inline
  std::size_t memory_consumption (const long double)
  {
    return sizeof(long double);
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
    std::size_t mem = sizeof(std::vector<T>);
    const unsigned int n = static_cast<unsigned int>(v.size());
    for (unsigned int i=0; i<n; ++i)
      mem += memory_consumption(v[i]);
    mem += (v.capacity() - n)*sizeof(T);
    return mem;
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



  inline
  std::size_t memory_consumption (const std::vector<int> &v)
  {
    return (v.capacity() * sizeof(int) +
            sizeof(v));
  }



  inline
  std::size_t memory_consumption (const std::vector<double> &v)
  {
    return (v.capacity() * sizeof(double) +
            sizeof(v));
  }



  inline
  std::size_t memory_consumption (const std::vector<float> &v)
  {
    return (v.capacity() * sizeof(float) +
            sizeof(v));
  }



  inline
  std::size_t memory_consumption (const std::vector<char> &v)
  {
    return (v.capacity() * sizeof(char) +
            sizeof(v));
  }



  inline
  std::size_t memory_consumption (const std::vector<unsigned char> &v)
  {
    return (v.capacity() * sizeof(unsigned char) +
            sizeof(v));
  }



  template <typename T>
  inline
  std::size_t memory_consumption (const std::vector<T *> &v)
  {
    return (v.capacity() * sizeof(T *) +
            sizeof(v));
  }



  template <typename A, typename B>
  inline
  std::size_t memory_consumption (const std::pair<A,B> &p)
  {
    return (memory_consumption(p.first) +
            memory_consumption(p.second));
  }



  template <typename T>
  inline
  std::size_t
  memory_consumption (const T *const)
  {
    return sizeof(T *);
  }



  template <typename T>
  inline
  std::size_t
  memory_consumption (T *const)
  {
    return sizeof(T *);
  }



  inline
  std::size_t
  memory_consumption (void *const)
  {
    return sizeof(void *);
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
  memory_consumption (const T &t)
  {
    return t.memory_consumption();
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
