//----------------------------  symmetric_tensor.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor.h  ---------------------------
#ifndef __deal2__symmetric_tensor_h
#define __deal2__symmetric_tensor_h


#include <base/tensor.h>


template <int rank, int dim> class SymmetricTensor;
template <int dim> class SymmetricTensor<2,dim>;


namespace internal
{
  namespace SymmetricTensor
  {
    namespace Rank2Accessors
    {

                                       /**
                                        * Switch type to select a tensor of
                                        * rank 2 and dimension <tt>dim</tt>,
                                        * switching on whether the tensor
                                        * should be constant or not.
                                        */
      template <int dim, bool constness>
      struct Types;

                                       /**
                                        * Switch type to select a tensor of
                                        * rank 2 and dimension <tt>dim</tt>,
                                        * switching on whether the tensor
                                        * should be constant or not.
                                        *
                                        * Specialization for constant tensors.
                                        */
      template <int dim>
      struct Types<dim,true>
      {
          typedef
          const typename ::SymmetricTensor<2,dim>::StorageType
          base_tensor_type;

          typedef double reference;
      };

                                       /**
                                        * Switch type to select a tensor of
                                        * rank 2 and dimension <tt>dim</tt>,
                                        * switching on whether the tensor
                                        * should be constant or not.
                                        *
                                        * Specialization for non-constant
                                        * tensors.
                                        */
      template <int dim>
      struct Types<dim,false>
      {
          typedef
          typename ::SymmetricTensor<2,dim>::StorageType
          base_tensor_type;

          typedef double &reference;
      };


                                       /**
                                        * Accessor class to access the
                                        * elements of individual rows in a
                                        * symmetric tensor. Since the elements
                                        * of symmetric tensors are not stored
                                        * as in a table, the accessors are a
                                        * little more involved.
                                        *
                                        * @author Wolfgang Bangerth, 2005
                                        */
      template <int dim, bool constness>
      class RowAccessor 
      {
        public:
                                           /**
                                            * Import which tensor we work on.
                                            */
          typedef
          typename Types<dim,constness>::base_tensor_type
          base_tensor_type;

                                           /**
                                            * The type of a reference to an
                                            * individual element of the
                                            * symmetric tensor. If the tensor
                                            * is constant, we can only return
                                            * a value instead of a reference.
                                            */
          typedef typename Types<dim,constness>::reference reference;

                                           /**
                                            * Constructor. Take the tensor to
                                            * access as well as the row we
                                            * point to as arguments.
                                            */
          RowAccessor (const base_tensor_type  &tensor,
                       const unsigned int  row);

                                           /**
                                            * Return a reference to an element
                                            * of this row (if we point to a
                                            * non-const tensor), or the value
                                            * of the element (in case this is
                                            * a constant tensor).
                                            */
          reference operator[] (const unsigned int column) const;
          
        private:
                                           /**
                                            * Reference to the tensor we
                                            * access.
                                            */
          const base_tensor_type &base_tensor;

                                           /**
                                            * Index of the row we access.
                                            */
          const unsigned int row;

                                           /**
                                            * Make the symmetric tensor
                                            * classes a friend, since they are
                                            * the only ones who can create
                                            * objects like this.
                                            */
          template <int,int> class ::SymmetricTensor;
      };

    }
  }
}

          


/**
 * Provide a class that stores symmetric tensors of rank 2 efficiently,
 * i.e. only store half of the off-diagonal elements of the full tensor.
 *
 * Using this tensor class for objects of rank 2 has advantages over
 * matrices in many cases since the dimension is known to the compiler
 * as well as the location of the data. It is therefore possible to
 * produce far more efficient code than for matrices with
 * runtime-dependent dimension.
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
class SymmetricTensor<2,dim>
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
    static const unsigned int rank      = 2;

                                     /**
                                      * Number of independent components of a
                                      * symmetric tensor of rank 2. We store
                                      * only the upper right half of it. This
                                      * information is probably of little
                                      * interest to all except the accessor
                                      * classes that need it.
                                      */
    static const unsigned int
    n_tensor_components = (dim*dim + dim)/2;

                                     /**
                                      * Declare the type in which we actually
                                      * store the data. This information is
                                      * probably of little interest to all
                                      * except the accessor classes that need
                                      * it. In particular, you shouldn't make
                                      * any assumptions about the storage
                                      * format in your application programs.
                                      */
    typedef Tensor<1,n_tensor_components> StorageType;
    
    
                                     /**
                                      * Default constructor. Creates a zero
                                      * tensor.
                                      */
    SymmetricTensor ();

                                     /**
                                      * Constructor. Generate a symmetric
                                      * tensor from a general one. Assumes
                                      * that @p t is already symmetric, but
                                      * this is not checked: we simply copy
                                      * only a subset of elements.
                                      */
    SymmetricTensor (const Tensor<2,dim> &t);

				     /**
				      *  Assignment operator.
				      */
    SymmetricTensor & operator = (const SymmetricTensor &);

				     /**
				      *  Test for equality of two tensors.
				      */
    bool operator == (const SymmetricTensor &) const;

    				     /**
				      *  Test for inequality of two tensors.
				      */
    bool operator != (const SymmetricTensor &) const;

				     /**
				      *  Add another tensor.
				      */
    SymmetricTensor & operator += (const SymmetricTensor &);
    
				     /**
				      *  Subtract another tensor.
				      */
    SymmetricTensor & operator -= (const SymmetricTensor &);

				     /**
				      *  Scale the tensor by <tt>factor</tt>,
				      *  i.e. multiply all components by
				      *  <tt>factor</tt>.
				      */
    SymmetricTensor & operator *= (const double factor);

				     /**
				      *  Scale the vector by
				      *  <tt>1/factor</tt>.
				      */
    SymmetricTensor & operator /= (const double factor);

				     /**
				      *  Add two tensors. If possible, you
				      *  should use <tt>operator +=</tt>
				      *  instead since this does not need the
				      *  creation of a temporary.
				      */
    SymmetricTensor   operator + (const SymmetricTensor &s) const;

				     /**
				      *  Subtract two tensors. If possible,
				      *  you should use <tt>operator -=</tt>
				      *  instead since this does not need the
				      *  creation of a temporary.
				      */
    SymmetricTensor   operator - (const SymmetricTensor &s) const;

				     /**
				      * Unary minus operator. Negate all
				      * entries of a tensor.
				      */
    SymmetricTensor   operator - () const;

                                     /**
                                      * Scalar product between two symmetric
                                      * tensors. It is the contraction
                                      * <tt>a<sub>ij</sub>b<sub>ij</sub></tt>
                                      * over all indices <tt>i,j</tt>. While
                                      * it is possible to define other scalar
                                      * products (and associated induced
                                      * norms), this one seems to be the most
                                      * appropriate one.
                                      */
    double operator * (const SymmetricTensor &s) const;
    
                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for constant tensors.
                                      */
    internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,true>
    operator [] (const unsigned int row) const;

                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for non-constant tensors.
                                      */
    internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,false>
    operator [] (const unsigned int row);
    
                                     /**
                                      * Return the Frobenius-norm of a tensor,
                                      * i.e. the square root of the sum of
                                      * squares of all entries. This norm is
                                      * induced by the scalar product defined
                                      * above for two symmetric tensors.
                                      */
    double norm () const;
    
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
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object.
				      */
    static unsigned int memory_consumption ();

    
  private:
                                     /**
                                      * Data storage for a symmetric tensor.
                                      */
    StorageType data;
};



// ------------------------- inline functions ------------------------

namespace internal
{
  namespace SymmetricTensor
  {
    namespace Rank2Accessors
    {
      template <int dim, bool constness>
      RowAccessor<dim,constness>::
      RowAccessor (const base_tensor_type &base_tensor,
                   const unsigned int row)
                      :
                      base_tensor (base_tensor),
                      row (row)
      {
        Assert (row < dim, ExcIndexRange (row, 0, dim));
      }



      template <int dim, bool constness>
      typename RowAccessor<dim,constness>::reference
      RowAccessor<dim,constness>::
      operator[] (const unsigned int column) const
      {
        Assert (column < dim, ExcIndexRange (column, 0, dim));

                                         // first treat the main diagonal
                                         // elements, which are stored
                                         // consecutively at the beginning
        if (row == column)
          return base_tensor[row];

                                         // the rest is messier and requires a
                                         // few switches. if someone has a
                                         // better idea, help is welcome
        switch (dim)
          {
            case 2:
                  Assert (((row==1) && (column==0)) || ((row==0) && (column==1)),
                          ExcInternalError());
                  return base_tensor[2];

            case 3:
                  if (((row==0) && (column==1)) ||
                      ((row==1) && (column==0)))
                    return base_tensor[3];
                  else if (((row==0) && (column==2)) ||
                           ((row==2) && (column==0)))
                    return base_tensor[4];
                  else if (((row==1) && (column==2)) ||
                           ((row==2) && (column==1)))
                    return base_tensor[5];
                  else
                    Assert (false, ExcInternalError());

            default:
                  Assert (false, ExcNotImplemented());
          }

        Assert (false, ExcInternalError());
        return 0;
      }
    }
  }
}



template <int dim>
inline
SymmetricTensor<2,dim>::SymmetricTensor ()
{}



template <>
inline
SymmetricTensor<2,2>::SymmetricTensor (const Tensor<2,2> &t)
{
  Assert (t[0][1] == t[1][0], ExcInternalError());

  data[0] = t[0][0];
  data[1] = t[1][1];
  data[2] = t[0][1];
}



template <>
inline
SymmetricTensor<2,3>::SymmetricTensor (const Tensor<2,3> &t)
{
  Assert (t[0][1] == t[1][0], ExcInternalError());
  Assert (t[0][2] == t[2][0], ExcInternalError());
  Assert (t[1][2] == t[2][1], ExcInternalError());
  
  data[0] = t[0][0];
  data[1] = t[1][1];
  data[2] = t[2][2];
  data[3] = t[0][1];
  data[4] = t[0][2];
  data[5] = t[1][2];
}


template <int dim>
inline
SymmetricTensor<2,dim> &
SymmetricTensor<2,dim>::operator = (const SymmetricTensor<2,dim> &t)
{
  data = t.data;
  return *this;
}



template <int dim>
inline
bool
SymmetricTensor<2,dim>::operator == (const SymmetricTensor<2,dim> &t) const
{
  return data == t.data;
}



template <int dim>
inline
bool
SymmetricTensor<2,dim>::operator != (const SymmetricTensor<2,dim> &t) const
{
  return data != t.data;
}



template <int dim>
inline
SymmetricTensor<2,dim> &
SymmetricTensor<2,dim>::operator += (const SymmetricTensor<2,dim> &t)
{
  data += t.data;
  return *this;
}



template <int dim>
inline
SymmetricTensor<2,dim> &
SymmetricTensor<2,dim>::operator -= (const SymmetricTensor<2,dim> &t)
{
  data -= t.data;
  return *this;
}



template <int dim>
inline
SymmetricTensor<2,dim> &
SymmetricTensor<2,dim>::operator *= (const double d)
{
  data *= d;
  return *this;
}



template <int dim>
inline
SymmetricTensor<2,dim> &
SymmetricTensor<2,dim>::operator /= (const double d)
{
  data /= d;
  return *this;
}



template <int dim>
inline
SymmetricTensor<2,dim>
SymmetricTensor<2,dim>::operator + (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data += t.data;
  return tmp;
}



template <int dim>
inline
SymmetricTensor<2,dim>
SymmetricTensor<2,dim>::operator - (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data -= t.data;
  return tmp;
}



template <int dim>
inline
SymmetricTensor<2,dim>
SymmetricTensor<2,dim>::operator - () const
{
  SymmetricTensor tmp = *this;
  tmp.data = -tmp.data;
  return tmp;
}



template <int dim>
inline
void
SymmetricTensor<2,dim>::clear ()
{
  data.clear ();
}



template <int dim>
inline
unsigned int
SymmetricTensor<2,dim>::memory_consumption ()
{
  return StorageType::memory_consumption ();
}



template <int dim>
double
SymmetricTensor<2,dim>::operator * (const SymmetricTensor &s) const
{
  double t = 0;
  unsigned int i=0;
  for (; i<dim; ++i)
    t += data[i] * s.data[i];

  for (; i<n_tensor_components; ++i)
    t += 2 * data[i] * s.data[i];

  return t;
}
  


template <int dim>
internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,true>
SymmetricTensor<2,dim>::operator [] (const unsigned int row) const
{
  return
    internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,true> (data, row);
}



template <int dim>
internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,false>
SymmetricTensor<2,dim>::operator [] (const unsigned int row)
{
  return
    internal::SymmetricTensor::Rank2Accessors::RowAccessor<dim,false> (data, row);
}



template <>
double
SymmetricTensor<2,1>::norm () const
{
  return std::sqrt(data[0]*data[0]);
}



template <>
double
SymmetricTensor<2,2>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + 2*data[2]*data[2]);
}



template <>
double
SymmetricTensor<2,3>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] +
                   2*data[3]*data[3] + 2*data[4]*data[4] + 2*data[5]*data[5]);
}


/* ----------------- Non-member functions operating on tensors. ------------ */

/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==2</tt>.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,2> &t)
{
  return (t[0][0] * t[1][1] - 2*t[0][1]*t[0][1]);
}




/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==3</tt>.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,3> &t)
{
				   // in analogy to general tensors, but
				   // there's something to be simplified for
				   // the present case
  return ( t[0][0]*t[1][1]*t[2][2]
	   -t[0][0]*t[1][2]*t[1][2]
	   -t[1][1]*t[0][2]*t[0][2]
	   -t[2][2]*t[0][1]*t[0][1]
	   +2*t[0][1]*t[0][2]*t[1][2] );
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of
 * its diagonal entries.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
double trace (const SymmetricTensor<2,dim> &d)
{
  double t=0;
  for (unsigned int i=0; i<dim; ++i)
    t += d[i][i];
  return t;
}



/**
 * Return the transpose of the given symmetric tensor. Since we are working
 * with symmetric objects, the transpose is of course the same as the original
 * tensor. This function mainly exists for compatibility with the Tensor
 * class.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<2,dim>
transpose (const SymmetricTensor<2,dim> &t)
{
  return t;
}


/**
 * Multiplication of a symmetric tensor of general rank with a scalar double
 * from the right.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator * (const SymmetricTensor<rank,dim> &t,
	    const double            factor)
{
  SymmetricTensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar double
 * from the left.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator * (const double            factor,
	    const SymmetricTensor<rank,dim> &t)
{
  SymmetricTensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Division of a symmetric tensor of general rank by a scalar double.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator / (const SymmetricTensor<rank,dim> &t,
	    const double            factor)
{
  SymmetricTensor<rank,dim> tt = t;
  tt /= factor;
  return tt;
}



#endif
