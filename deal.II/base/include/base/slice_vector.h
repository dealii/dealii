//---------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------
#ifndef __deal2__slice_vector_h
#define __deal2__slice_vector_h

#include <base/config.h>
#include <base/exceptions.h>

/**
 * Filter a range out of any object having a random access operator
 * <tt>[] (unsigned int)</tt>.
 *
 * The use of this object is straight forward. It reduplicates the
 * random access operator of the <tt>VECTOR</tt> and adds an offset to
 * every index.
 *
 * Some precautions have to be taken if it is used for a constant
 * vector: the SliceVector object has to be constant, too. The
 * appropriate initalization sequence is like this:
 *
 * @code
 *   void f(const std::vector<int>& v)
 *   {
 *     const SliceVector slice(v,...);
 *     ...
 *   }
 * @endcode
 *
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
class SliceVector
{
  public:
    SliceVector(VECTOR& v);
    SliceVector(VECTOR& v,
		unsigned int start,
		unsigned int length);

    unsigned int size() const;

    typename VECTOR::reference operator[] (unsigned int);
    typename VECTOR::const_reference operator[] (unsigned int) const;
  private:
    VECTOR& v;
    unsigned int start;
    unsigned int length;
};


/**
 * Helper function for creating temporary objects without typing
 * template arguments.
 *
 * @relates SliceVector
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
inline
const SliceVector<VECTOR>
make_slice (VECTOR& v)
{
  const SliceVector<VECTOR> r(v);
  return r;
}



/**
 * Helper function for creating temporary objects without typing
 * template arguments.
 *
 * @relates SliceVector
 * @author Guido Kanschat, 2004
 */
template <class VECTOR>
inline
const SliceVector<VECTOR>
make_slice (VECTOR& v, unsigned int start, unsigned int length)
{
  const SliceVector<VECTOR> r(v, start, length);
  return r;
}




//-------------- member functions --------------------//

template <class VECTOR>
inline
SliceVector<VECTOR>::SliceVector(VECTOR& v)
		: v(v), start(0), length(v.size())
{}


template <class VECTOR>
inline
SliceVector<VECTOR>::SliceVector(VECTOR& v,
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
SliceVector<VECTOR>::size() const
{
  return length;
}


template <class VECTOR>
inline
typename VECTOR::reference
SliceVector<VECTOR>::operator[](unsigned int i)
{
  Assert ((i<length), ExcIndexRange(i, 0, length));
  
  return v[start+i];
}


template <class VECTOR>
inline
typename VECTOR::const_reference
SliceVector<VECTOR>::operator[](unsigned int i) const
{
  Assert ((i<length), ExcIndexRange(i, 0, length));
  
  return v[start+i];
}


#endif
