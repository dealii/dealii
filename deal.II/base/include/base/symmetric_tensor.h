//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__symmetric_tensor_h
#define __deal2__symmetric_tensor_h


#include <base/tensor.h>
#include <base/table_indices.h>


template <int rank, int dim> class SymmetricTensor;


namespace internal
{
                                   /**
                                    * A namespace for classes that are
                                    * internal to how the SymmetricTensor
                                    * class works.
                                    */
  namespace SymmetricTensor
  {
                                     /**
                                      * Declaration of typedefs for the type
                                      * of data structures which are used to
                                      * store symmetric tensors. For example,
                                      * for rank-2 symmetric tensors, we use a
                                      * flat vector to store all the
                                      * elements. On the other hand, symmetric
                                      * rank-4 tensors are mappings from
                                      * symmetric rank-2 tensors into
                                      * symmetric rank-2 tensors, so they can
                                      * be represented as matrices, etc.
                                      *
                                      * This information is probably of little
                                      * interest to all except the accessor
                                      * classes that need it. In particular,
                                      * you shouldn't make any assumptions
                                      * about the storage format in your
                                      * application programs.
                                      */
    template <int rank, int dim>
    struct StorageType;

                                     /**
                                      * Specialization of StorageType for
                                      * rank-2 tensors.
                                      */
    template <int dim>
    struct StorageType<2,dim> 
    {
                                         /**
                                          * Number of independent components of a
                                          * symmetric tensor of rank 2. We store
                                          * only the upper right half of it.
                                          */
        static const unsigned int
        n_independent_components = (dim*dim + dim)/2;

                                         /**
                                          * Declare the type in which we actually
                                          * store the data.
                                          */
        typedef Tensor<1,n_independent_components> base_tensor_type;
    };



                                     /**
                                      * Specialization of StorageType for
                                      * rank-4 tensors.
                                      */
    template <int dim>
    struct StorageType<4,dim> 
    {
                                         /**
                                          * Number of independent components
                                          * of a symmetric tensor of rank
                                          * 2. Since rank-4 tensors are
                                          * mappings between such objects, we
                                          * need this information.
                                          */
        static const unsigned int
        n_rank2_components = (dim*dim + dim)/2;

                                         /**
                                          * Declare the type in which we
                                          * actually store the data. Symmetric
                                          * rank-4 tensors are mappings
                                          * between symmetric rank-2 tensors,
                                          * so we can represent the data as a
                                          * matrix if we represent the rank-2
                                          * tensors as vectors.
                                          */
        typedef Tensor<2,n_rank2_components> base_tensor_type;
    };
    
    

                                     /**
                                      * Switch type to select a tensor of
                                      * rank 2 and dimension <tt>dim</tt>,
                                      * switching on whether the tensor
                                      * should be constant or not.
                                      */
    template <int rank, int dim, bool constness>
    struct AccessorTypes;

                                     /**
                                      * Switch type to select a tensor of
                                      * rank 2 and dimension <tt>dim</tt>,
                                      * switching on whether the tensor
                                      * should be constant or not.
                                      *
                                      * Specialization for constant tensors.
                                      */
    template <int rank, int dim>
    struct AccessorTypes<rank, dim,true>
    {
        typedef const ::SymmetricTensor<rank,dim> tensor_type;

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
    template <int rank, int dim>
    struct AccessorTypes<rank,dim,false>
    {
        typedef ::SymmetricTensor<rank,dim> tensor_type;

        typedef double &reference;
    };


    template <int rank, int dim, bool constness>
    class Accessor;
    
                                     /**
                                      * Accessor class to access the elements
                                      * of individual rows in a symmetric
                                      * tensor of rank 2. Since the elements
                                      * of symmetric tensors are not stored as
                                      * in a table, the accessors are a little
                                      * more involved. However, for tensors of
                                      * rank 2 they are still relatively
                                      * simple in that an accessor is created
                                      * by the SymmetricTensor class with the
                                      * first access to <tt>operator[]</tt>;
                                      * the accessor thereby points to a row
                                      * of the tensor. Calling
                                      * <tt>operator[]</tt> on the accessor
                                      * then selects an entry of this
                                      * row. Note that if this entry is not
                                      * actually stored, then the transpose
                                      * entry is chosen as that is guaranteed
                                      * to be stored.
                                      *
                                      * @author Wolfgang Bangerth, 2005
                                      */
    template <int dim, bool constness>
    class Accessor<2,dim,constness> 
    {
      public:
                                         /**
                                          * Import which tensor we work on.
                                          */
        typedef
        typename AccessorTypes<2,dim,constness>::tensor_type
        tensor_type;

                                         /**
                                          * The type of a reference to an
                                          * individual element of the
                                          * symmetric tensor. If the tensor
                                          * is constant, we can only return
                                          * a value instead of a reference.
                                          */
        typedef typename AccessorTypes<2,dim,constness>::reference reference;

                                         /**
                                          * Constructor. Take the tensor to
                                          * access as well as the row we
                                          * point to as arguments.
                                          */
        Accessor (tensor_type        &tensor,
                  const unsigned int  row);

                                         /**
                                          * Return a reference to an element
                                          * of this row (if we point to a
                                          * non-const tensor), or the value
                                          * of the element (in case this is
                                          * a constant tensor).
                                          */
        reference operator[] (const unsigned int column);
          
      private:
                                         /**
                                          * Reference to the tensor we
                                          * access.
                                          */
        tensor_type &tensor;

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

          


/**
 * Provide a class that stores symmetric tensors of rank 2,4,... efficiently,
 * i.e. only store those off-diagonal elements of the full tensor that are not
 * redundant. For example, for symmetric 2x2 tensors, this would be the
 * elements 11, 22, and 12, while the element 21 is equal to the 12 element.
 *
 * Using this class for symmetric tensors of rank 2 has advantages over
 * matrices in many cases since the dimension is known to the compiler as well
 * as the location of the data. It is therefore possible to produce far more
 * efficient code than for matrices with runtime-dependent dimension. It is
 * also more efficient than using the more general <tt>Tensor</tt> class,
 * since less elements are stored, and the class automatically makes sure that
 * the tensor represents a symmetric object.
 *
 * For tensors of higher rank, the savings in storage are even higher. For
 * example for the 3x3x3x3 tensors of rank 4, only 36 instead of the full 81
 * entries have to be stored.
 *
 * Tensors of rank 4 are considered symmetric if they are operators mapping
 * symmetric rank-2 tensors onto symmetric rank-2 tensors. This entails
 * certain symmetry properties on the elements in their 4-dimensional index
 * space.
 *
 * Symmetric tensors are most often used in structural and fluid mechanics,
 * where strains and stresses are usually symmetric tensors, and the
 * stress-strain relationship is given by a symmetric rank-4 tensor.
 *
 * Note that symmetric tensors only exist with even numbers of indices. In
 * other words, the only objects that you can use are
 * <tt>SymmetricTensor<2,dim></tt>, <tt>SymmetricTensor<4,dim></tt>, etc, but
 * <tt>SymmetricTensor<1,dim></tt> and <tt>SymmetricTensor<3,dim></tt> do not
 * exist and their use will most likely lead to compiler errors.
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int rank, int dim>
class SymmetricTensor
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
                                      * Return a read-write reference
                                      * to the indicated element.
                                      */
    double & operator() (const TableIndices<rank> &indices);
  
                                     /**
                                      * Return the value of the
                                      * indicated element as a
                                      * read-only reference.
                                      *
                                      * We return the requested value
                                      * as a constant reference rather
                                      * than by value since this
                                      * object may hold data types
                                      * that may be large, and we
                                      * don't know here whether
                                      * copying is expensive or not.
                                      */
    double operator() (const TableIndices<rank> &indices) const;

                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for constant tensors.
                                      */
    internal::SymmetricTensor::Accessor<rank,dim,true>
    operator [] (const unsigned int row) const;

                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for non-constant tensors.
                                      */
    internal::SymmetricTensor::Accessor<rank,dim,false>
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
    typename internal::SymmetricTensor::StorageType<2,dim>::base_tensor_type data;
};



// ------------------------- inline functions ------------------------

namespace internal
{
  namespace SymmetricTensor
  {
    template <int dim, bool constness>
    Accessor<2,dim,constness>::
    Accessor (tensor_type       &tensor,
              const unsigned int  row)
                    :
                    tensor (tensor),
                    row (row)
    {}



    template <int dim, bool constness>
    typename Accessor<2,dim,constness>::reference
    Accessor<2,dim,constness>::
    operator[] (const unsigned int column)
    {
      return tensor(TableIndices<2> (row, column));
    }
  }
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>::SymmetricTensor ()
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


template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator = (const SymmetricTensor<rank,dim> &t)
{
  data = t.data;
  return *this;
}



template <int rank, int dim>
inline
bool
SymmetricTensor<rank,dim>::operator == (const SymmetricTensor<rank,dim> &t) const
{
  return data == t.data;
}



template <int rank, int dim>
inline
bool
SymmetricTensor<rank,dim>::operator != (const SymmetricTensor<rank,dim> &t) const
{
  return data != t.data;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator += (const SymmetricTensor<rank,dim> &t)
{
  data += t.data;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator -= (const SymmetricTensor<rank,dim> &t)
{
  data -= t.data;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator *= (const double d)
{
  data *= d;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator /= (const double d)
{
  data /= d;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator + (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data += t.data;
  return tmp;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator - (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data -= t.data;
  return tmp;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator - () const
{
  SymmetricTensor tmp = *this;
  tmp.data = -tmp.data;
  return tmp;
}



template <int rank, int dim>
inline
void
SymmetricTensor<rank,dim>::clear ()
{
  data.clear ();
}



template <int rank, int dim>
inline
unsigned int
SymmetricTensor<rank,dim>::memory_consumption ()
{
  return
    internal::SymmetricTensor::StorageType<rank,dim>::memory_consumption ();
}



template <>
double
SymmetricTensor<2,1>::operator * (const SymmetricTensor<2,1> &s) const
{
  return data[0] * s.data[0];
}



template <>
double
SymmetricTensor<2,2>::operator * (const SymmetricTensor<2,2> &s) const
{
  return (data[0] * s.data[0] +
          data[1] * s.data[1] +
          2*data[2] * s.data[2]);
}



template <>
double
SymmetricTensor<2,3>::operator * (const SymmetricTensor<2,3> &s) const
{
  return (data[0] * s.data[0] +
          data[1] * s.data[1] +
          data[2] * s.data[2] +
          2*data[3] * s.data[3] +
          2*data[4] * s.data[4] +
          2*data[5] * s.data[5]);
}



template <>
double &
SymmetricTensor<2,1>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0];
}



template <>
double
SymmetricTensor<2,1>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0];
}



template <>
double &
SymmetricTensor<2,2>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches. at least for the 2x2 case it
                                   // is reasonably simple
  Assert (((indices[0]==1) && (indices[1]==0)) ||
          ((indices[0]==0) && (indices[1]==1)),
          ExcInternalError());  
  return data[2];
}



template <>
double
SymmetricTensor<2,2>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches. at least for the 2x2 case it
                                   // is reasonably simple
  Assert (((indices[0]==1) && (indices[1]==0)) ||
          ((indices[0]==0) && (indices[1]==1)),
          ExcInternalError());  
  return data[2];
}



template <>
double &
SymmetricTensor<2,3>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches, but simpler if we just sort
                                   // our indices
  TableIndices<2> sorted_indices (indices);
  sorted_indices.sort ();
  
  if ((sorted_indices[0]==0) && (sorted_indices[1]==1))
    return data[3];
  else if ((sorted_indices[0]==0) && (sorted_indices[1]==2))
    return data[4];
  else if ((sorted_indices[0]==1) && (sorted_indices[1]==2))
    return data[5];
  else
    Assert (false, ExcInternalError());

  static double dummy_but_referenceable = 0;
  return dummy_but_referenceable;
}



template <>
double
SymmetricTensor<2,3>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches, but simpler if we just sort
                                   // our indices
  TableIndices<2> sorted_indices (indices);
  sorted_indices.sort ();
  
  if ((sorted_indices[0]==0) && (sorted_indices[1]==1))
    return data[3];
  else if ((sorted_indices[0]==0) && (sorted_indices[1]==2))
    return data[4];
  else if ((sorted_indices[0]==1) && (sorted_indices[1]==2))
    return data[5];
  else
    Assert (false, ExcInternalError());

  static double dummy_but_referenceable = 0;
  return dummy_but_referenceable;
}



template <int rank, int dim>
internal::SymmetricTensor::Accessor<rank,dim,true>
SymmetricTensor<rank,dim>::operator [] (const unsigned int row) const
{
  return
    internal::SymmetricTensor::Accessor<rank,dim,true> (*this, row);
}



template <int rank, int dim>
internal::SymmetricTensor::Accessor<rank,dim,false>
SymmetricTensor<rank,dim>::operator [] (const unsigned int row)
{
  return
    internal::SymmetricTensor::Accessor<rank,dim,false> (*this, row);
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
template <int rank, int dim>
double trace (const SymmetricTensor<rank,dim> &d)
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
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
transpose (const SymmetricTensor<rank,dim> &t)
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
