// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __dealii__vector_view_h
#define __dealii__vector_view_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector.h>

#include <cstdio>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Vectors
 *@{
 */

/**
 * View of a numerical vector of data. This class provides an
 * interface compatible with the Vector<double> class (from which it
 * is inherited), that allows fast access to locations of memory
 * already allocated with arrays of type Number.
 *
 * This is in the same style of the vector view in the Trilinos
 * library.
 *
 * You should consider using the VectorView object ONLY when ALL of
 * the following requirements are met:
 *
 * 1. Your application requires a Vector<Number> object.
 *
 * 2. All you have at your disposal is a Number* pointer.
 *
 * 3. You are ABSOLUTELY SURE that the above pointer points to a
 * valid area of memory of the correct size.
 *
 * 4. You really believe that making a copy of that memory would be
 * too expensive.
 *
 * 5. You really know what you are doing.
 *
 * Notice that NO CHECKS are performed on the actual memory, and if
 * you try to access illegal areas of memory, your computer will
 * suffer from it. Once again, use this class ONLY if you know exactly
 * what you are doing.
 *
 * Two constructors are provided. One for read-write access, and one
 * for read only access, and you are allowed to use this class on
 * objects of type const Number*. However you should be aware of the
 * fact that the constness of the array pointed to is ignored, which
 * means that you should only use the const constructor when the
 * actual object you are constructing is itself a constant object.
 * As a corollary, you will be allowed to call even functions of
 * the base class that change data of the array; this being a violation
 * of the C++ type model, you should make sure that this only happens
 * if it is really valid and, in general, if
 * you know what you are doing.
 *
 * Since this class does not own the memory that you are accessing,
 * you have to make sure that the lifetime of the section of memory you
 * are viewing is longer than this object. No attempt is made to
 * ensure that this is the case.
 *
 * An example usage of this class is the following:
 *
 @code
 // Create an array of length 5;
 double * array = new double[5];
 // Now create a view of the above array that is compatible with the
 // Vector<double> class
 VectorView<double> view(5, array);

 view(1) = 4;

 // The following line should output 4.
 cout << array[1] << endl;

 // If debug mode is on, then the following triggers an execption:
 view(6) = 4;

 // But notice that no checks are performed, so this is legal but WILL
 // NOT work
 VectorView<double> wrong_view(10, array);

 // Now no assert will be thrown if you type wrong_view(6), but most
 // likely a seg fault will occur.
 view(6) = 4;

 // Notice that this construction is legal. It will create a copy of
 // the array.
 const Vector<double> const_copy(view);

 // Now this is the correct way to instantiate a constant view of the
 // above vector:
 const VectorView<double> correct_const_copy_view(5, const_copy.begin());

 // While this will compile, BUT WILL NOT COMPLAIN if you try to write
 // on it!
 VectorView<double> wrong_const_copy_view(5, const_copy.begin());

 // Now writing to elements of wrong_const_copy_view is allowed, and
 // will change the same memory as the const_copy object.
 wrong_const_copy_view(1) = 5;

 if(copy_view(1) == wrong_const_copy_view(1)) cout << "Tautology";

 @endcode
 *
 *
 * @note Instantiations for this template are provided for
 * <tt>@<float@>, @<double@>, @<long double@>,
 * @<std::complex@<float@>@>, @<std::complex@<double@>@>,
 * @<std::complex@<long double@>@></tt>; others can be generated in
 * application programs (see the section on @ref Instantiations in the
 * manual).
 *
 * @author Luca Heltai, 2009
 */
template<typename Number>
class VectorView : public Vector<Number>
{
public:

  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Read write constructor. Takes the size
   * of the vector, just like the standard
   * one, but the data is picked starting
   * from the location of the pointer @p
   * ptr.
   */
  VectorView(const size_type new_size, Number *ptr);

  /**
   * The constant constructor is the same
   * as above, however you will not be able
   * to access the data for write access.
   *
   * You should only use this class by
   * constructing it as a const
   * VectorView<double>(size, ptr) object.
   *
   * Undefined behavior will occur if you
   * construct it as a non const object or
   * attempt to write on it.
   */
  VectorView(const size_type new_size, const Number *ptr);

  /**
   * This desctructor will only reset the
   * internal sizes and the internal
   * pointers, but it will NOT clear the
   * memory. */
  ~VectorView();

  /**
   * The reinit function of this object has
   * a behavior which is different from the
   * one of the base class. VectorView does
   * not handle memory, and you should not
   * attempt to resize the memory that is
   * pointed to by this object using the
   * reinit function. You can, however,
   * resize the view that you have of the
   * original object. Notice that it is
   * your own responsibility to ensure that
   * the memory you are pointing to is big
   * enough.
   *
   * Similarly to what happens in the base
   * class, if fast is false, then the
   * entire content of the vector is set to
   * 0, otherwise the content of the memory
   * is left unchanged.
   *
   * Notice that the following snippet of
   * code may not produce what you expect:
   *
   * @code
   * // Create a vector of length 1.
   * Vector<double> long_vector(1);
   *
   * // Make a view of it
   * VectorView<double> view_of_long_vector(1, long_vector.begin());
   *
   * // Resize the original vector to a bigger size
   * long_vector.reinit(100);
   *
   * // And the view, leaving the memory untouched
   * view_of_long_vector.reinit(100, true);
   * @endcode
   *
   * In the above case, the
   * Vector<double>::reinit method is
   * called, and a NEW area of memory is
   * reserved, possibly not starting at the
   * same place as before. Hoever, the
   * VectorView<double> object keeps
   * pointing to the same old area. After
   * the two reinits, any call to
   * view_of_long_vector(i), with i>0 might
   * cause an attempt to access invalid
   * areas of memory, or might function
   * properly, depending on whether or not
   * the system was able to allocate some
   * memory consecutively after the
   * original allocation.
   *
   * In any case, you should not rely on
   * this behavior, and you should only
   * call this reinit function if you
   * really know what you are doing.
   */
  virtual void reinit (const size_type N,
                       const bool         fast=false);

  /** This reinit function is
      equivalent to constructing a
      new object with the given
      size, starting from the
      pointer ptr. */
  void reinit(const size_type N, Number *ptr);

  /** This reinit function is
      equivalent to constructing a
      new object with the given
      size, starting from the
      pointer ptr. The same
      considerations made for the
      constructor apply here. */
  void reinit(const size_type N, const Number *ptr);

  /**
   * This function is here to prevent
   * memory corruption. It should never be
   * called, and will throw an exception if
   * you try to do so. */
  virtual void swap (Vector<Number> &v);
};



/*@}*/
/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

template<typename Number>
inline
VectorView<Number>::VectorView(const size_type new_size, Number *ptr)
{
  this->vec_size      = new_size;
  this->max_vec_size  = new_size;
  this->val           = ptr;
}



template<typename Number>
inline
VectorView<Number>::VectorView(const size_type new_size, const Number *ptr)
{
  this->vec_size      = new_size;
  this->max_vec_size  = new_size;
  this->val           = const_cast<Number *>(ptr);
}



template<typename Number>
inline
VectorView<Number>::~VectorView()
{
  // avoid that the base class releases
  // memory it doesn't own
  this->vec_size = 0;
  this->max_vec_size = 0;
  this->val = 0;
}


template<typename Number>
inline
void VectorView<Number>::reinit(const size_type N, const bool fast)
{
  this->vec_size = N;
  this->max_vec_size = N;
  if (fast == false)
    Vector<Number>::operator=(static_cast<Number>(0));
}


template<typename Number>
inline
void VectorView<Number>::reinit(const size_type new_size, Number *ptr)
{
  this->vec_size      = new_size;
  this->max_vec_size  = new_size;
  this->val           = ptr;
}


template<typename Number>
inline
void VectorView<Number>::reinit(const size_type new_size, const Number *ptr)
{
  this->vec_size      = new_size;
  this->max_vec_size  = new_size;
  this->val           = const_cast<Number *>(ptr);
}


template<typename Number>
inline
void VectorView<Number>::swap(Vector<Number> &)
{
  AssertThrow(false, ExcMessage("Cant' swap a VectorView with a Vector!"));
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
