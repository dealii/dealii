//----------------------------  tensor_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_base.h  ---------------------------
#ifndef __deal2__tensor_base_h
#define __deal2__tensor_base_h


// single this file out from tensor.h, since we want to derive Point<dim>
// from Tensor<1,dim>. However, the point class will not need all the
// tensor stuff, so we don't want the whole tensor package to be included
// everytime we use a point.


#include <base/config.h>
#include <base/exceptions.h>
#include <vector>

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif

template <typename number> class Vector;
template <int dim> class Point;

// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank, int dim> class Tensor;
template <int dim> class Tensor<1,dim>;


/**
 * This class is a specialized version of the <tt>Tensor<rank,dim></tt> class.
 * It handles tensors with one index, i.e. vectors, of fixed dimension and
 * provides the basis for the functionality needed for tensors of higher rank.
 *
 * Within deal.II, the distinction between this class and its derived class
 * <tt>Point</tt> is that we use the <tt>Point</tt> class mainly to denote the
 * points that make up geometric objects. As such, they have a small number of
 * additional operations over general tensors of rank 1 for which we use the
 * <tt>Tensor<1,dim></tt> class. In particular, there is a distance() function
 * to compute the Euclidian distance between two points in space.
 *
 * However, the <tt>Point</tt> class is really only used where the coordinates
 * of an object can be thought to possess the dimension of a length. For all
 * other uses, such as the gradient of a scalar function (which is a tensor of
 * rank 1, or vector, with as many elements as a point object, but with
 * different physical units), we use the <tt>Tensor<1,dim></tt> class.
 */
template <int dim>
class Tensor<1,dim>
{
  public:
				     /**
				      * Provide a way to get the
				      * dimension of an object without
				      * explicit knowledge of it's
				      * data type. Implementation is
				      * this way instead of providing
				      * a function <tt>dimension()</tt>
				      * because now it is possible to
				      * get the dimension at compile
				      * time without the expansion and
				      * preevaluation of an inlined
				      * function; the compiler may
				      * therefore produce more
				      * efficient code and you may use
				      * this value to declare other
				      * data types.
				      */
    static const unsigned int dimension = dim;

				     /**
				      * Publish the rank of this tensor to
				      * the outside world.
				      */
    static const unsigned int rank      = 1;

				     /**
				      * Type of stored objects. This
				      * is a double for a rank 1 tensor.
				      */

    typedef double value_type;
    
				     /**
				      * Declare an array type which can
				      * be used to initialize statically
				      * an object of this type.
				      *
				      * Avoid warning about zero-sized
				      * array for <tt>dim==0</tt> by
				      * choosing lunatic value that is
				      * likely to overflow memory
				      * limits.
				      */
    typedef double array_type[(dim!=0) ? dim : 100000000];
    
				     /**
				      * Constructor. Initialize all entries
				      * to zero if <tt>initialize==true</tt>; this
				      * is the default behaviour.
				      */
    explicit Tensor (const bool initialize = true);

				     /**
				      * Copy constructor, where the
				      * data is copied from a C-style
				      * array.
				      */
    Tensor (const array_type &initializer);
    
    				     /**
				      * Copy constructor.
				      */
    Tensor (const Tensor<1,dim> &);

				     /**
				      * Read access to the <tt>index</tt>th
				      * coordinate.
				      *
				      * Note that the derived
				      * <tt>Point</tt> class also provides
				      * access through the <tt>()</tt>
				      * operator for
				      * backcompatibility.
				      */
    double   operator [] (const unsigned int index) const;

    				     /**
				      * Read and write access to the
				      * <tt>index</tt>th coordinate.
				      *
				      * Note that the derived
				      * <tt>Point</tt> class also provides
				      * access through the <tt>()</tt>
				      * operator for
				      * backcompatibility.
				      */
    double & operator [] (const unsigned int index);

				     /**
				      * Assignment operator.
				      */
    Tensor<1,dim> & operator = (const Tensor<1,dim> &);

				     /**
				      * Test for equality of two
				      * tensors.
				      */
    bool operator == (const Tensor<1,dim> &) const;

    				     /**
				      * Test for inequality of two
				      * tensors.
				      */
    bool operator != (const Tensor<1,dim> &) const;

				     /**
				      * Add another vector, i.e. move
				      * this point by the given
				      * offset.
				      */
    Tensor<1,dim> & operator += (const Tensor<1,dim> &);
    
				     /**
				      * Subtract another vector.
				      */
    Tensor<1,dim> & operator -= (const Tensor<1,dim> &);

				     /**
				      * Scale the vector by
				      * <tt>factor</tt>, i.e. multiply all
				      * coordinates by <tt>factor</tt>.
				      */
    Tensor<1,dim> & operator *= (const double factor);

				     /**
				      * Scale the vector by <tt>1/factor</tt>.
				      */
    Tensor<1,dim> & operator /= (const double factor);

				     /**
				      * Returns the scalar product of
				      * two vectors.
				      */
    double          operator * (const Tensor<1,dim> &) const;

				     /**
				      * Add two tensors. If possible,
				      * use <tt>operator +=</tt> instead
				      * since this does not need to
				      * copy a point at least once.
				      */
    Tensor<1,dim>   operator + (const Tensor<1,dim> &) const;

				     /**
				      * Subtract two tensors. If
				      * possible, use <tt>operator +=</tt>
				      * instead since this does not
				      * need to copy a point at least
				      * once.
				      */
    Tensor<1,dim>   operator - (const Tensor<1,dim> &) const;

				     /**
				      * Tensor with inverted entries.
				      */
    Tensor<1,dim>   operator - () const;
    
				     /**
				      * Reset all values to zero.
				      *
				      * Note that this is partly inconsistent
				      * with the semantics of the @p clear()
				      * member functions of the STL and of
				      * several other classes within deal.II
				      * which not only reset the values of
				      * stored elements to zero, but release
				      * all memory and return the object into
				      * a virginial state. However, since the
				      * size of objects of the present type is
				      * determined by its template parameters,
				      * resizing is not an option, and indeed
				      * the state where all elements have a
				      * zero value is the state right after
				      * construction of such an object.
				      */
    void clear ();

				     /**
				      * Fill a vector with all tensor elements.
				      *
				      * This function unrolls all
				      * tensor entries into a single,
				      * linearly numbered vector. As
				      * usual in C++, the rightmost
				      * index marches fastest.
				      */
    void unroll (Vector<double> &result) const;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object.
				      */
    static unsigned int memory_consumption ();

    				     /** @addtogroup Exceptions
				      * @{ */
    
                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcDimTooSmall,
                    int,
                    << "Given dimensions must be >= 1, but was " << arg1);

				     //@}
  private:
				     /**
				      * Store the values in a simple
				      * array.  For <tt>dim==0</tt> store
				      * one element, because otherways
				      * the compiler would choke.  We
				      * catch this case in the
				      * constructor to disallow the
				      * creation of such an object.
				      */
    double values[(dim!=0) ? (dim) : 1];

#ifdef DEAL_II_TEMPLATE_SPEC_ACCESS_WORKAROUND
  public:
#endif
				     /**
				      * Help function for unroll. If
				      * we have detected an access
				      * control bug in the compiler,
				      * this function is declared
				      * public, otherwise private. Do
				      * not attempt to use this
				      * function from outside in any
				      * case, even if it should be
				      * public for your compiler.
				      */
    void unroll_recursion (Vector<double> &result,
			   unsigned int   &start_index) const;
    
  private:
				     /**
				      * Make the following classes
				      * friends to this class. In
				      * principle, it would suffice if
				      * otherrank==2, but that is not
				      * possible in C++ at present.
				      *
				      * Also, it would be sufficient
				      * to make the function
				      * unroll_loops a friend, but
				      * that seems to be impossible as
				      * well.
				      */
    template <int otherrank, int otherdim>  friend class Tensor;

				     /**
				      * Point is allowed access to
				      * the coordinates. This is
				      * supposed to improve speed.
				      */
    friend class Point<dim>;
};


				 /**
				  *  Prints the values of this point in the
				  *  form <tt>x1 x2 x3 etc</tt>.
				  */
template <int dim>
std::ostream & operator << (std::ostream &out, const Tensor<1,dim> &p);

/// @if NoDoc

/*------------------------------- Inline functions: Tensor ---------------------------*/


template <int dim>
inline
Tensor<1,dim>::Tensor (const bool initialize)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  if (initialize)
    for (unsigned int i=0; i!=dim; ++i)
      values[i] = 0;
}



template <int dim>
inline
Tensor<1,dim>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}



template <int dim>
inline
Tensor<1,dim>::Tensor (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
}



template <>
inline
Tensor<1,0>::Tensor (const Tensor<1,0> &)
{
				   // at some places in the library,
				   // we have Point<0> for formal
				   // reasons (e.g., we sometimes have
				   // Quadrature<dim-1> for faces, so
				   // we have Quadrature<0> for dim=1,
				   // and then we have Point<0>). To
				   // avoid warnings in the above
				   // function that the loop end check
				   // always fails, we implement this
				   // function here
}



template <int dim>
inline
double Tensor<1,dim>::operator [] (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <int dim>
inline
double & Tensor<1,dim>::operator [] (const unsigned int index)
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <>
inline
Tensor<1,0> & Tensor<1,0>::operator = (const Tensor<1,0> &)
{
				   // at some places in the library,
				   // we have Point<0> for formal
				   // reasons (e.g., we sometimes have
				   // Quadrature<dim-1> for faces, so
				   // we have Quadrature<0> for dim=1,
				   // and then we have Point<0>). To
				   // avoid warnings in the above
				   // function that the loop end check
				   // always fails, we implement this
				   // function here
  return *this;
}



template <>
inline
Tensor<1,1> & Tensor<1,1>::operator = (const Tensor<1,1> &p)
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  values[0] = p.values[0];
  return *this;
}



template <>
inline
Tensor<1,2> & Tensor<1,2>::operator = (const Tensor<1,2> &p)
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  values[0] = p.values[0];
  values[1] = p.values[1];
  return *this;
}



template <>
inline
Tensor<1,3> & Tensor<1,3>::operator = (const Tensor<1,3> &p)
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  values[0] = p.values[0];
  values[1] = p.values[1];
  values[2] = p.values[2];
  return *this;
}



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator = (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
  return *this;
}



template <int dim>
inline
bool Tensor<1,dim>::operator == (const Tensor<1,dim> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i])
      return false;
  return true;
}



template <int dim>
inline
bool Tensor<1,dim>::operator != (const Tensor<1,dim> &p) const
{
  return !((*this) == p);
}



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator += (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
}



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator -= (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
}



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator *= (const double s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
}



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator /= (const double s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
}



template <>
inline
double Tensor<1,1>::operator * (const Tensor<1,1> &p) const
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  return (values[0] * p.values[0]);
}



template <>
inline
double Tensor<1,2>::operator * (const Tensor<1,2> &p) const
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  return (values[0] * p.values[0] +
	  values[1] * p.values[1]);
}



template <>
inline
double Tensor<1,3>::operator * (const Tensor<1,3> &p) const
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  return (values[0] * p.values[0] +
	  values[1] * p.values[1] +
	  values[2] * p.values[2]);
}



template <int dim>
inline
double Tensor<1,dim>::operator * (const Tensor<1,dim> &p) const
{
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += values[i] * p.values[i];
  return q;
}



template <int dim>
inline
Tensor<1,dim> Tensor<1,dim>::operator + (const Tensor<1,dim> &p) const
{
  return (Tensor<1,dim>(*this) += p);
}



template <int dim>
inline
Tensor<1,dim> Tensor<1,dim>::operator - (const Tensor<1,dim> &p) const
{
  return (Tensor<1,dim>(*this) -= p);
}



template <int dim>
inline
Tensor<1,dim> Tensor<1,dim>::operator - () const
{
  Tensor<1,dim> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -values[i];
  return result;
}



template <int dim>
inline
void Tensor<1,dim>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
}



template <int dim>
inline
unsigned int
Tensor<1,dim>::memory_consumption ()
{
  return sizeof(Tensor<1,dim>);
}



/**
 * Output operator for tensors of rank 1. Print the elements
 * consecutively, with a space in between.
 */
template <int dim>
inline
std::ostream & operator << (std::ostream &out, const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}



/**
 * Output operator for tensors of rank 1 and dimension 1. This is
 * implemented specialized from the general template in order to avoid
 * a compiler warning that the loop is empty.
 */
inline
std::ostream & operator << (std::ostream &out, const Tensor<1,1> &p)
{
  out << p[0];

  return out;
}



/**
 * Multiplication of a tensor of rank 1 with a scalar double from the right.
 */
template <int dim>
inline
Tensor<1,dim>
operator * (const Tensor<1,dim> &t,
	    const double         factor)
{
  Tensor<1,dim> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



/**
 * Multiplication of a tensor of rank 1 with a scalar double from the left.
 */
template <int dim>
inline
Tensor<1,dim>
operator * (const double         factor,
	    const Tensor<1,dim> &t)
{
  Tensor<1,dim> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



/**
 * Division of a tensor of rank 1 by a scalar double.
 */
template <int dim>
inline
Tensor<1,dim>
operator / (const Tensor<1,dim> &t,
	    const double         factor)
{
  Tensor<1,dim> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] / factor;
  return tt;
}

/// @endif

#endif


