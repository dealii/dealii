//----------------------------  tensor_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
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

#include <base/exceptions.h>
#include <iostream>
#include <vector>

template <typename number> class Vector;


// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank, int dim> class Tensor;


/**
 * This class is a specialized version of the #Tensor<rank,dim># class.
 * It handles tensors with one index, i.e. vectors, of fixed dimension
 * and offers the functionality needed for tensors of higher rank.
 *
 * In many cases, you will want to use the more specialized #Point# class
 * which acts as a tensor of rank one but has more functionality.
 */
template <int dim>
class Tensor<1,dim> {
      public:
				     /**
				      * Provide a way to get the
				      * dimension of an object without
				      * explicit knowledge of it's
				      * data type. Implementation is
				      * this way instead of providing
				      * a function #dimension()#
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

#if (__GNUC__==2) && (__GNUC_MINOR__ < 95)  
				     /**
				      * Unnecessary variable, only used in
				      * the declaration of #array_type#.
				      * This variable should be omitted,
				      * but egcs1.1.2 chokes on that so
				      * we need it presently.
				      */
    static const unsigned int array_size = ((dim!=0) ? (dim) : 1);

				     /**
				      * Declare an array type which can
				      * be used to initialize statically
				      * an object of this type.
				      *
				      * Use the same condition for #dim==0#
				      * as in the #TensorBase<1,dim># class.
				      */
    typedef double array_type[array_size];
#else

				     /**
				      * Declare an array type which can
				      * be used to initialize statically
				      * an object of this type.
				      *
				      * Use the same condition for #dim==0#
				      * as in the #TensorBase<1,dim># class.
				      */
    typedef double array_type[(dim!=0) ? dim : 1];
#endif
    
				     /**
				      * Constructor. Initialize all entries
				      * to zero if #initialize==true#; this
				      * is the default behaviour.
				      */
    explicit Tensor (const bool initialize = true);

				     /**
				      * Copy constructor, where the data is
				      * copied from a C-style array.
				      */
    Tensor (const array_type &initializer);
    
    				     /**
				      *  Copy constructor.
				      */
    Tensor (const Tensor<1,dim> &);

				     /**
				      *  Read access to the #index#th coordinate.
				      *
				      * Note that the derived #Point# class also
				      * provides access through the #()#
				      * operator for backcompatibility.
				      */
    double   operator [] (const unsigned int index) const;

    				     /**
				      *  Read and write access to the #index#th
				      *  coordinate.
				      *
				      * Note that the derived #Point# class also
				      * provides access through the #()#
				      * operator for backcompatibility.
				      */
    double & operator [] (const unsigned int index);

				     /**
				      *  Assignment operator.
				      */
    Tensor<1,dim> & operator = (const Tensor<1,dim> &);

				     /**
				      *  Test for equality of two points.
				      */
    bool operator == (const Tensor<1,dim> &) const;

    				     /**
				      *  Test for inequality of two points.
				      */
    bool operator != (const Tensor<1,dim> &) const;

				     /**
				      *  Add another vector, i.e. move this point by
				      *  the given offset.
				      */
    Tensor<1,dim> & operator += (const Tensor<1,dim> &);
				     /**
				      *  Subtract another vector.
				      */
    Tensor<1,dim> & operator -= (const Tensor<1,dim> &);

				     /**
				      *  Scale the vector by #factor#, i.e. multiply
				      *  all coordinates by #factor#.
				      */
    Tensor<1,dim> & operator *= (const double &factor);

				     /**
				      *  Scale the vector by #1/factor#.
				      */
    Tensor<1,dim> & operator /= (const double &factor);

				     /**
				      *  Returns the scalar product of two vectors.
				      */
    double          operator * (const Tensor<1,dim> &) const;

				     /**
				      *  Add two tensors. If possible, use
				      *  #operator +=# instead since this does not
				      *  need to copy a point at least once.
				      */
    Tensor<1,dim>   operator + (const Tensor<1,dim> &) const;

				     /**
				      *  Subtract two tensors. If possible, use
				      *  #operator +=# instead since this does not
				      *  need to copy a point at least once.
				      */
    Tensor<1,dim>   operator - (const Tensor<1,dim> &) const;

				     /**
				      * Reset all values to zero.
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

  protected:
				     /**
				      * Store the values in a simple array.
				      * For #dim==0# store one element, because
				      * otherways the compiler would choke.
				      * We catch this case in the constructor
				      * to disallow the creation of such
				      * an object.
				      */
    double values[(dim!=0) ? (dim) : 1];

				     /**
				      * Help function for unroll.
				      */
    void unroll_recursion (Vector<double> &result,
			   unsigned int   &start_index) const;

				     // make the following classes
				     // friends to this class. in principle,
				     // it would suffice if otherrank==2,
				     // but then the compiler complains
				     // that this be an explicit specialization
				     // which is not what we want
				     //
				     // also, it would be sufficient to make
				     // the function unroll_loops a friend,
				     // but that seems to be impossible as well.
    template<int otherrank, int otherdim>  friend class Tensor;
};

				 /**
				  * Exception
				  */
DeclException2(ExcWrongVectorSize, int, int, << "Tensor has " << arg1
	       << " entries, but vector has size " << arg2);
DeclException1 (ExcDimTooSmall,
		int,
		<< "Given dimensions must be >= 1, but was " << arg1);


				 /**
				  *  Prints the values of this point in the
				  *  form #x1 x2 x3 etc#.
				  */
template <int dim>
ostream & operator << (ostream &out, const Tensor<1,dim> &p);


/*------------------------------- Inline functions: Tensor ---------------------------*/
// in the following, declarations of operators look rather weird because they are
// split over several lines. this is necessary, however, because otherwise doc++
// gets really confused about this all... :-(


template <int dim>
inline
Tensor<1,dim>::Tensor (const bool initialize)
{
  Assert (dim>0, ExcDimTooSmall(dim));

  if (initialize)
    for (unsigned int i=0; i!=dim; ++i)
      values[i] = 0;
};


template <int dim>
inline
Tensor<1,dim>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = initializer[i];
};


template <int dim>
inline
Tensor<1,dim>::Tensor (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
};


template <int dim>
inline
double Tensor<1,dim>::operator [] (const unsigned int index) const
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
};


template <int dim>
inline
double & Tensor<1,dim>::operator [] (const unsigned int index)
{
  Assert (index<dim, ExcIndexRange (index, 0, dim));
  return values[index];
};


template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator = (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
  return *this;
};


template <int dim>
inline
bool Tensor<1,dim>::operator == (const Tensor<1,dim> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i]) return false;
  return true;
};


template <int dim>
inline
bool Tensor<1,dim>::operator != (const Tensor<1,dim> &p) const
{
  return !((*this) == p);
};


template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator += (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
};


template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator -= (const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
};


template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator *= (const double &s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
};


template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator /= (const double &s)
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
};


template <int dim>
inline
double Tensor<1,dim>::operator * (const Tensor<1,dim> &p) const
{
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += values[i] * p.values[i];
  return q;
};


template <int dim>
inline
Tensor<1,dim> Tensor<1,dim>::operator + (const Tensor<1,dim> &p) const
{
  return (Tensor<1,dim>(*this) += p);
};


template <int dim>
inline
Tensor<1,dim> Tensor<1,dim>::operator - (const Tensor<1,dim> &p) const
{
  return (Tensor<1,dim>(*this) -= p);
};


template <int dim>
inline
void Tensor<1,dim>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
};


template <int dim>
inline
ostream & operator << (ostream &out, const Tensor<1,dim> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
};


template <>
inline
ostream & operator << (ostream &out, const Tensor<1,1> &p)
{
  out << p[0];

  return out;
};


#endif
