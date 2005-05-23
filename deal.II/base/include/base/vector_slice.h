//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__vector_slice_h
#define __deal2__vector_slice_h

#include <base/config.h>
#include <base/exceptions.h>

/**
 * Filter a range out of any object having a random access operator
 * <tt>[] (unsigned int)</tt> and a function <tt>size() const</tt>.
 *
 * The use of this object is straight forward. It reduplicates the
 * random access operator of the <tt>VECTOR</tt> and adds an offset to
 * every index.
 *
 * Some precautions have to be taken if it is used for a constant
 * vector: the VectorSlice object has to be constant, too. The
 * appropriate initalization sequence is like this:
 *
 * @code
 *   void f(const std::vector<int>& v)
 *   {
 *     const VectorSlice slice(v,...);
 *     ...
 *   }
 * @endcode
 *
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
class VectorSlice
{
  public:
				     /**
				      * Construct a vector slice
				      * containing the whole
				      * vector. Comes handy, if you
				      * did not want to have a slice
				      * at all, but the function you
				      * call wants it: just put in the
				      * vector itself as argument and
				      * let this constructor make a
				      * slice for you.
				      */
    VectorSlice(VECTOR& v);
				     /**
				      * The real constructor for a
				      * vector slice, allowing you to
				      * specify the start index and
				      * the length of the slice.
				      */
    VectorSlice(VECTOR& v,
		unsigned int start,
		unsigned int length);

				     /**
				      * Return the length of the slice
				      * using the same interface as
				      * <tt>std::vector</tt>.
				      */
    unsigned int size() const;

				     /**
				      * Access an element of the slice
				      * using the same interface as
				      * <tt>std::vector</tt>.
				      */
    typename VECTOR::reference operator[] (unsigned int);
    
				     /**
				      * Access an element of a
				      * constant slice using the same
				      * interface as
				      * <tt>std::vector</tt>.
				      */
    typename VECTOR::const_reference operator[] (unsigned int) const;
  private:
				     /**
				      * The vector we extract from.
				      */
    VECTOR& v;
				     /**
				      * The start index of the slice.
				      */
    const unsigned int start;
				     /**
				      * The length of the slice.
				      */
    const unsigned int length;
};


/**
 * Helper function for creating temporary objects without typing
 * template arguments.
 *
 * @relates VectorSlice
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
inline
const VectorSlice<const VECTOR>
make_slice (VECTOR& v)
{
  const VectorSlice<const VECTOR> r(v);
  return r;
}



/**
 * Helper function for creating temporary objects without typing
 * template arguments.
 *
 * @relates VectorSlice
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
inline
const VectorSlice<const VECTOR>
make_slice (VECTOR& v, unsigned int start, unsigned int length)
{
  const VectorSlice<const VECTOR> r(v, start, length);
  return r;
}




//---------------------------------------------------------------------------

template <class VECTOR>
inline
VectorSlice<VECTOR>::VectorSlice(VECTOR& v)
		: v(v), start(0), length(v.size())
{}


template <class VECTOR>
inline
VectorSlice<VECTOR>::VectorSlice(VECTOR& v,
			 unsigned int start,
			 unsigned int length)
		: v(v), start(start), length(length)
{
  Assert((start+length<=v.size()),
	 ExcIndexRange(length, 0, v.size()-start+1));
}


template <class VECTOR>
inline
unsigned int
VectorSlice<VECTOR>::size() const
{
  return length;
}


template <class VECTOR>
inline
typename VECTOR::reference
VectorSlice<VECTOR>::operator[](unsigned int i)
{
  Assert ((i<length), ExcIndexRange(i, 0, length));
  
  return v[start+i];
}


template <class VECTOR>
inline
typename VECTOR::const_reference
VectorSlice<VECTOR>::operator[](unsigned int i) const
{
  Assert ((i<length), ExcIndexRange(i, 0, length));
  
  return v[start+i];
}


#endif
