//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tensor_base_h
#define __deal2__tensor_base_h


// single this file out from tensor.h, since we want to derive
// Point<dim,Number> from Tensor<1,dim,Number>. However, the point class will
// not need all the tensor stuff, so we don't want the whole tensor package to
// be included everytime we use a point.


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <vector>

#include <cmath>

DEAL_II_NAMESPACE_OPEN
// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif

template <typename number> class Vector;

// forward declare Point and Tensor. This is the first definition of these
// classes and here we set the default Number type to double (this means that
// this file must be included when using something like Tensor<1,dim>, and
// Point and Tensor must not be forward declared without the number type
// specified)
template <int dim, typename Number=double> class Point;

// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank_, int dim, typename Number=double> class Tensor;
template <int dim, typename Number> class Tensor<0,dim,Number>;
template <int dim, typename Number> class Tensor<1,dim,Number>;



/**
 * This class is a specialized version of the <tt>Tensor<rank,dim,Number></tt>
 * class. It handles tensors of rank zero, i.e. scalars. The second template
 * argument is ignored.
 *
 * This class exists because in some cases we want to construct
 * objects of type Tensor@<spacedim-dim,dim,Number@>, which should expand to
 * scalars, vectors, matrices, etc, depending on the values of the
 * template arguments @p dim and @p spacedim. We therefore need a
 * class that acts as a scalar (i.e. @p Number) for all purposes but
 * is part of the Tensor template family.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 2009
 */
template <int dim, typename Number>
class Tensor<0,dim,Number>
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
    static const unsigned int rank      = 0;

				     /**
				      * Type of stored objects. This
				      * is a Number for a rank 1 tensor.
				      */

    typedef Number value_type;

				     /**
				      * Declare a type that has holds
				      * real-valued numbers with the same
				      * precision as the template argument to
				      * this class. For std::complex<number>,
				      * this corresponds to type number, and
				      * it is equal to Number for all other
				      * cases. See also the respective field
				      * in Vector<Number>.
				      *
				      * This typedef is used to
				      * represent the return type of
				      * norms.
				      */
    typedef typename numbers::NumberTraits<Number>::real_type real_type;

				     /**
				      * Constructor. Set to zero.
				      */
    Tensor ();

				     /**
				      * Copy constructor, where the
				      * data is copied from a C-style
				      * array.
				      */
    Tensor (const value_type &initializer);

    				     /**
				      * Copy constructor.
				      */
    Tensor (const Tensor<0,dim,Number> &);

				     /**
				      * Conversion to Number. Since
				      * rank-0 tensors are scalars,
				      * this is a natural operation.
				      */
    operator Number () const;

				     /**
				      * Conversion to Number. Since
				      * rank-0 tensors are scalars,
				      * this is a natural operation.
				      *
				      * This is the non-const
				      * conversion operator that
				      * returns a writable reference.
				      */
    operator Number& ();

				     /**
				      * Assignment operator.
				      */
    Tensor<0,dim,Number> & operator = (const Tensor<0,dim,Number> &);

    				     /**
				      * Assignment operator.
				      */
    Tensor<0,dim,Number> & operator = (const Number d);

				     /**
				      * Test for equality of two
				      * tensors.
				      */
    bool operator == (const Tensor<0,dim,Number> &) const;

    				     /**
				      * Test for inequality of two
				      * tensors.
				      */
    bool operator != (const Tensor<0,dim,Number> &) const;

				     /**
				      * Add another vector, i.e. move
				      * this point by the given
				      * offset.
				      */
    Tensor<0,dim,Number> & operator += (const Tensor<0,dim,Number> &);

				     /**
				      * Subtract another vector.
				      */
    Tensor<0,dim,Number> & operator -= (const Tensor<0,dim,Number> &);

				     /**
				      * Scale the vector by
				      * <tt>factor</tt>, i.e. multiply all
				      * coordinates by <tt>factor</tt>.
				      */
    Tensor<0,dim,Number> & operator *= (const Number factor);

				     /**
				      * Scale the vector by <tt>1/factor</tt>.
				      */
    Tensor<0,dim,Number> & operator /= (const Number factor);

				     /**
				      * Returns the scalar product of
				      * two vectors.
				      */
    Number                 operator * (const Tensor<0,dim,Number> &) const;

				     /**
				      * Add two tensors. If possible,
				      * use <tt>operator +=</tt> instead
				      * since this does not need to
				      * copy a point at least once.
				      */
    Tensor<0,dim,Number>   operator + (const Tensor<0,dim,Number> &) const;

				     /**
				      * Subtract two tensors. If
				      * possible, use <tt>operator +=</tt>
				      * instead since this does not
				      * need to copy a point at least
				      * once.
				      */
    Tensor<0,dim,Number>   operator - (const Tensor<0,dim,Number> &) const;

				     /**
				      * Tensor with inverted entries.
				      */
    Tensor<0,dim,Number>   operator - () const;

                                     /**
                                      * Return the Frobenius-norm of a
                                      * tensor, i.e. the square root
                                      * of the sum of squares of all
                                      * entries. For the present case
                                      * of rank-1 tensors, this equals
                                      * the usual
                                      * <tt>l<sub>2</sub></tt> norm of
                                      * the vector.
                                      */
    real_type norm () const;

                                     /**
                                      * Return the square of the
                                      * Frobenius-norm of a tensor,
                                      * i.e. the square root of the
                                      * sum of squares of all entries.
				      *
				      * This function mainly exists
				      * because it makes computing the
				      * norm simpler recursively, but
				      * may also be useful in other
				      * contexts.
                                      */
    real_type norm_square () const;

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
                                      * Only tensors with a positive
				      * dimension are implemented. This
				      * exception is thrown by the
				      * constructor if the template
				      * argument <tt>dim</tt> is zero or
				      * less.
				      *
				      * @ingroup Exceptions
                                      */
    DeclException1 (ExcDimTooSmall,
                    int,
                    << "dim must be positive, but was " << arg1);

                     /**
                      * Read or write the data of this object to or
                      * from a stream for the purpose of serialization
                      */
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version);

  private:
				     /**
				      * The value of this scalar object.
				      */
    Number value;
};



/**
 * This class is a specialized version of the <tt>Tensor<rank,dim,Number></tt>
 * class.  It handles tensors with one index, i.e. vectors, of fixed dimension
 * and provides the basis for the functionality needed for tensors of higher
 * rank.
 *
 * Within deal.II, the distinction between this class and its derived class
 * <tt>Point</tt> is that we use the <tt>Point</tt> class mainly to denote the
 * points that make up geometric objects. As such, they have a small number of
 * additional operations over general tensors of rank 1 for which we use the
 * <tt>Tensor<1,dim,Number></tt> class. In particular, there is a distance()
 * function to compute the Euclidian distance between two points in space.
 *
 * However, the <tt>Point</tt> class is really only used where the coordinates
 * of an object can be thought to possess the dimension of a length. For all
 * other uses, such as the gradient of a scalar function (which is a tensor of
 * rank 1, or vector, with as many elements as a point object, but with
 * different physical units), we use the <tt>Tensor<1,dim,Number></tt> class.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 1998-2005
 */
template <int dim,typename Number>
class Tensor<1,dim,Number>
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
				      * is a Number for a rank 1 tensor.
				      */

    typedef Number value_type;

				     /**
				      * Declare a type that has holds
				      * real-valued numbers with the same
				      * precision as the template argument to
				      * this class. For std::complex<number>,
				      * this corresponds to type number, and
				      * it is equal to Number for all other
				      * cases. See also the respective field
				      * in Vector<Number>.
				      *
				      * This typedef is used to
				      * represent the return type of
				      * norms.
				      */
    typedef typename numbers::NumberTraits<Number>::real_type real_type;

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
    typedef Number array_type[(dim!=0) ? dim : 100000000];

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
    Tensor (const Tensor<1,dim,Number> &);

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
    Number   operator [] (const unsigned int index) const;

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
    Number & operator [] (const unsigned int index);

				     /**
				      * Assignment operator.
				      */
    Tensor<1,dim,Number> & operator = (const Tensor<1,dim,Number> &);

    				     /**
				      * This operator assigns a scalar
				      * to a tensor. To avoid
				      * confusion with what exactly it
				      * means to assign a scalar value
				      * to a tensor, zero is the only
				      * value allowed for <tt>d</tt>,
				      * allowing the intuitive
				      * notation <tt>t=0</tt> to reset
				      * all elements of the tensor to
				      * zero.
				      */
    Tensor<1,dim,Number> & operator = (const Number d);

				     /**
				      * Test for equality of two
				      * tensors.
				      */
    bool operator == (const Tensor<1,dim,Number> &) const;

    				     /**
				      * Test for inequality of two
				      * tensors.
				      */
    bool operator != (const Tensor<1,dim,Number> &) const;

				     /**
				      * Add another vector, i.e. move
				      * this point by the given
				      * offset.
				      */
    Tensor<1,dim,Number> & operator += (const Tensor<1,dim,Number> &);

				     /**
				      * Subtract another vector.
				      */
    Tensor<1,dim,Number> & operator -= (const Tensor<1,dim,Number> &);

				     /**
				      * Scale the vector by
				      * <tt>factor</tt>, i.e. multiply all
				      * coordinates by <tt>factor</tt>.
				      */
    Tensor<1,dim,Number> & operator *= (const Number factor);

				     /**
				      * Scale the vector by <tt>1/factor</tt>.
				      */
    Tensor<1,dim,Number> & operator /= (const Number factor);

				     /**
				      * Returns the scalar product of
				      * two vectors.
				      */
    Number                 operator * (const Tensor<1,dim,Number> &) const;

				     /**
				      * Add two tensors. If possible,
				      * use <tt>operator +=</tt> instead
				      * since this does not need to
				      * copy a point at least once.
				      */
    Tensor<1,dim,Number>   operator + (const Tensor<1,dim,Number> &) const;

				     /**
				      * Subtract two tensors. If
				      * possible, use <tt>operator +=</tt>
				      * instead since this does not
				      * need to copy a point at least
				      * once.
				      */
    Tensor<1,dim,Number>   operator - (const Tensor<1,dim,Number> &) const;

				     /**
				      * Tensor with inverted entries.
				      */
    Tensor<1,dim,Number>   operator - () const;

                                     /**
                                      * Return the Frobenius-norm of a
                                      * tensor, i.e. the square root
                                      * of the sum of squares of all
                                      * entries. For the present case
                                      * of rank-1 tensors, this equals
                                      * the usual
                                      * <tt>l<sub>2</sub></tt> norm of
                                      * the vector.
                                      */
    real_type norm () const;

                                     /**
                                      * Return the square of the
                                      * Frobenius-norm of a tensor,
                                      * i.e. the square root of the
                                      * sum of squares of all entries.
				      *
				      * This function mainly exists
				      * because it makes computing the
				      * norm simpler recursively, but
				      * may also be useful in other
				      * contexts.
                                      */
    real_type norm_square () const;

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
    template <typename Number2>
    void unroll (Vector<Number2> &result) const;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object.
				      */
    static std::size_t memory_consumption ();

                                     /**
                                      * Only tensors with a positive
				      * dimension are implemented. This
				      * exception is thrown by the
				      * constructor if the template
				      * argument <tt>dim</tt> is zero or
				      * less.
				      *
				      * @ingroup Exceptions
                                      */
    DeclException1 (ExcDimTooSmall,
                    int,
                    << "dim must be positive, but was " << arg1);

                     /**
                      * Read or write the data of this object to or
                      * from a stream for the purpose of serialization
                      */
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version);

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
    Number values[(dim!=0) ? (dim) : 1];

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
    template <typename Number2>
    void unroll_recursion (Vector<Number2> &result,
			   unsigned int    &start_index) const;

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
    template <int otherrank, int otherdim, typename OtherNumber>  friend class dealii::Tensor;

				     /**
				      * Point is allowed access to
				      * the coordinates. This is
				      * supposed to improve speed.
				      */
    friend class Point<dim>;
};


				 /**
				  *  Prints the value of this scalar.
				  */
template <int dim,typename Number>
std::ostream & operator << (std::ostream &out, const Tensor<0,dim,Number> &p);

				 /**
				  *  Prints the values of this tensor in the
				  *  form <tt>x1 x2 x3 etc</tt>.
				  */
template <int dim,typename Number>
std::ostream & operator << (std::ostream &out, const Tensor<1,dim,Number> &p);


#ifndef DOXYGEN

/*---------------------------- Inline functions: Tensor<0,dim> ------------------------*/

template <int dim,typename Number>
inline
Tensor<0,dim,Number>::Tensor ()
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = 0;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::Tensor (const value_type &initializer)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = initializer;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::Tensor (const Tensor<0,dim,Number> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  value = p.value;
}




template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator Number () const
{
  return value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number>::operator Number & ()
{
  return value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator = (const Tensor<0,dim,Number> &p)
{
  value = p.value;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator = (const Number d)
{
  value = d;
  return *this;
}



template <int dim, typename Number>
inline
bool Tensor<0,dim,Number>::operator == (const Tensor<0,dim,Number> &p) const
{
  return (value == p.value);
}



template <int dim, typename Number>
inline
bool Tensor<0,dim,Number>::operator != (const Tensor<0,dim,Number> &p) const
{
  return !((*this) == p);
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator += (const Tensor<0,dim,Number> &p)
{
  value += p.value;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator -= (const Tensor<0,dim,Number> &p)
{
  value -= p.value;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator *= (const Number s)
{
  value *= s;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> & Tensor<0,dim,Number>::operator /= (const Number s)
{
  value /= s;
  return *this;
}



template <int dim, typename Number>
inline
Number Tensor<0,dim,Number>::operator * (const Tensor<0,dim,Number> &p) const
{
  return value*p.value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> Tensor<0,dim,Number>::operator + (const Tensor<0,dim,Number> &p) const
{
  return value+p.value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> Tensor<0,dim,Number>::operator - (const Tensor<0,dim,Number> &p) const
{
  return value-p.value;
}



template <int dim, typename Number>
inline
Tensor<0,dim,Number> Tensor<0,dim,Number>::operator - () const
{
  return -value;
}



template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
Tensor<0,dim,Number>::norm () const
{
  return numbers::NumberTraits<Number>::abs (value);
}



template <int dim, typename Number>
inline
typename Tensor<0,dim,Number>::real_type
Tensor<0,dim,Number>::norm_square () const
{
  return numbers::NumberTraits<Number>::abs_square (value);
}



template <int dim, typename Number>
inline
void Tensor<0,dim,Number>::clear ()
{
  value = 0;
}



template <int dim, typename Number>
template <class Archive>
inline
void Tensor<0,dim,Number>::serialize(Archive & ar, const unsigned int)
{
  ar & value;
}

/*---------------------------- Inline functions: Tensor<1,dim,Number> ------------------------*/


template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const bool initialize)
{
  if (initialize)
    for (unsigned int i=0; i!=dim; ++i)
      values[i] = 0;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const array_type &initializer)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number>::Tensor (const Tensor<1,dim,Number> &p)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
}



template <>
inline
Tensor<1,0,double>::Tensor (const Tensor<1,0,double> &)
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



template <int dim, typename Number>
inline
Number Tensor<1,dim,Number>::operator [] (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <int dim, typename Number>
inline
Number & Tensor<1,dim,Number>::operator [] (const unsigned int index)
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
}



template <>
inline
Tensor<1,0,double> & Tensor<1,0,double>::operator = (const Tensor<1,0,double> &)
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



template <int dim, typename Number>
inline
Tensor<1,dim,Number> &
Tensor<1,dim,Number>::operator = (const Tensor<1,dim,Number> &p)
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  switch (dim)
    {
    case 1:
      values[0] = p.values[0];
      break;
    case 2:
      values[0] = p.values[0];
      values[1] = p.values[1];
      break;
    case 3:
      values[0] = p.values[0];
      values[1] = p.values[1];
      values[2] = p.values[2];
      break;
    default:
      for (unsigned int i=0; i<dim; ++i)
	values[i] = p.values[i];
    }
  return *this;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> & Tensor<1,dim,Number>::operator = (const Number d)
{
  Assert (d==Number(0), ExcMessage ("Only assignment with zero is allowed"));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;

  return *this;
}



template <int dim, typename Number>
inline
bool Tensor<1,dim,Number>::operator == (const Tensor<1,dim,Number> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i])
      return false;
  return true;
}



template <>
inline
bool Tensor<1,0,double>::operator == (const Tensor<1,0,double> &) const
{
  return true;
}



template <int dim, typename Number>
inline
bool Tensor<1,dim,Number>::operator != (const Tensor<1,dim,Number> &p) const
{
  return !((*this) == p);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> & Tensor<1,dim,Number>::operator += (const Tensor<1,dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> & Tensor<1,dim,Number>::operator -= (const Tensor<1,dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> & Tensor<1,dim,Number>::operator *= (const Number s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> & Tensor<1,dim,Number>::operator /= (const Number s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
}



template <int dim, typename Number>
inline
Number
Tensor<1,dim,Number>::operator * (const Tensor<1,dim,Number> &p) const
{
				   // unroll by hand since this is a
				   // frequently called function and
				   // some compilers don't want to
				   // always unroll the loop in the
				   // general template
  switch (dim)
    {
    case 1:
      return (values[0] * p.values[0]);
      break;
    case 2:
      return (values[0] * p.values[0] +
	      values[1] * p.values[1]);
      break;
    case 3:
      return (values[0] * p.values[0] +
	      values[1] * p.values[1] +
	      values[2] * p.values[2]);
      break;
    default:
      Number q=0;
      for (unsigned int i=0; i<dim; ++i)
	q += values[i] * p.values[i];
      return q;
    }
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> Tensor<1,dim,Number>::operator + (const Tensor<1,dim,Number> &p) const
{
  return (Tensor<1,dim,Number>(*this) += p);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> Tensor<1,dim,Number>::operator - (const Tensor<1,dim,Number> &p) const
{
  return (Tensor<1,dim,Number>(*this) -= p);
}



template <int dim, typename Number>
inline
Tensor<1,dim,Number> Tensor<1,dim,Number>::operator - () const
{
  Tensor<1,dim,Number> result;
  for (unsigned int i=0; i<dim; ++i)
    result.values[i] = -values[i];
  return result;
}



template <int dim, typename Number>
inline
typename Tensor<1,dim,Number>::real_type
Tensor<1,dim,Number>::norm () const
{
  return std::sqrt (norm_square());
}



template <int dim, typename Number>
inline
typename Tensor<1,dim,Number>::real_type
Tensor<1,dim,Number>::norm_square () const
{
  real_type s = 0;
  for (unsigned int i=0; i<dim; ++i)
    s += numbers::NumberTraits<Number>::abs_square(values[i]);

  return s;
}



template <int dim, typename Number>
template <typename Number2>
inline
void
Tensor<1,dim,Number>::unroll (Vector<Number2> &result) const
{
  Assert (result.size()==dim,
          ExcDimensionMismatch(dim, result.size()));

  unsigned index = 0;
  unroll_recursion (result,index);
}



template<int dim, typename Number>
template <typename Number2>
inline
void
Tensor<1,dim,Number>::unroll_recursion (Vector<Number2> &result,
					unsigned int    &index) const
{
  for (unsigned i=0; i<dim; ++i)
    result(index++) = operator[](i);
}



template <int dim, typename Number>
inline
void Tensor<1,dim,Number>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
}



template <int dim, typename Number>
inline
std::size_t
Tensor<1,dim,Number>::memory_consumption ()
{
  return sizeof(Tensor<1,dim,Number>);
}


template <int dim, typename Number>
template <class Archive>
inline
void Tensor<1,dim,Number>::serialize(Archive & ar, const unsigned int)
{
  ar & values;
}
#endif // DOXYGEN


/**
 * Output operator for tensors of rank 0. Since such tensors are
 * scalars, we simply print this one value.
 *
 * @relates Tensor<0,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream & operator << (std::ostream &out, const Tensor<0,dim,Number> &p)
{
  out << static_cast<Number>(p);
  return out;
}



/**
 * Output operator for tensors of rank 1. Print the elements
 * consecutively, with a space in between.
 *
 * @relates Tensor<1,dim,Number>
 */
template <int dim, typename Number>
inline
std::ostream & operator << (std::ostream &out, const Tensor<1,dim,Number> &p)
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
 *
 * @relates Tensor<1,dim,Number>
 */
inline
std::ostream & operator << (std::ostream &out, const Tensor<1,1,double> &p)
{
  out << p[0];

  return out;
}



/**
 * Multiplication of a tensor of rank 1 with a scalar Number from the right.
 *
 * @relates Tensor<1,dim,Number>
 */
template <int dim, typename Number>
inline
Tensor<1,dim,Number>
operator * (const Tensor<1,dim,Number> &t,
	    const Number                factor)
{
  Tensor<1,dim,Number> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



/**
 * Multiplication of a tensor of rank 1 with a scalar Number from the left.
 *
 * @relates Tensor<1,dim,Number>
 */
template <int dim, typename Number>
inline
Tensor<1,dim,Number>
operator * (const Number                factor,
	    const Tensor<1,dim,Number> &t)
{
  Tensor<1,dim,Number> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



/**
 * Division of a tensor of rank 1 by a scalar Number.
 *
 * @relates Tensor<1,dim,Number>
 */
template <int dim, typename Number>
inline
Tensor<1,dim,Number>
operator / (const Tensor<1,dim,Number> &t,
	    const Number                factor)
{
  Tensor<1,dim,Number> tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] / factor;
  return tt;
}



/**
 * Multiplication of a tensor of rank 1 with a scalar double from the right.
 *
 * @relates Tensor<1,dim,Number>
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
 *
 * @relates Tensor<1,dim,Number>
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
 *
 * @relates Tensor<1,dim,Number>
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


DEAL_II_NAMESPACE_CLOSE

#endif


