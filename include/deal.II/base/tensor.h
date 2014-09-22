// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__tensor_h
#define __deal2__tensor_h


#include <deal.II/base/config.h>
#include <deal.II/base/tensor_base.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN
template <int rank_, int dim, typename Number> class Tensor;
template <int dim, typename Number> class Tensor<1,dim,Number>;

/**
 * Provide a general tensor class with an arbitrary rank, i.e. with
 * an arbitrary number of indices. The Tensor class provides an
 * indexing operator and a bit of infrastructure, but most
 * functionality is recursively handed down to tensors of rank 1 or
 * put into external templated functions, e.g. the <tt>contract</tt> family.
 *
 * Using this tensor class for objects of rank 2 has advantages over
 * matrices in many cases since the dimension is known to the compiler
 * as well as the location of the data. It is therefore possible to
 * produce far more efficient code than for matrices with
 * runtime-dependent dimension.
 *
 * This class provides an optional template argument for the type of the
 * underlying data. It defaults to @p double values. It can be used to base
 * tensors on @p float or @p complex numbers or any other data type that
 * implements basic arithmetic operations.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 1998-2005
 */
template <int rank_, int dim, typename Number>
class Tensor
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
  static const unsigned int rank      = rank_;

  /**
   * Number of independent components of a
   * tensor of current rank. This is dim times the
   * number of independent components of each sub-tensor.
   */
  static const unsigned int
  n_independent_components = Tensor<rank_-1,dim>::n_independent_components *dim;

  /**
   * Type of stored objects. This
   * is a tensor of lower rank.
   */
  typedef Tensor<rank_-1,dim,Number> value_type;

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
   * Declare an array type which
   * can be used to initialize an
   * object of this type
   * statically.
   */
  typedef typename Tensor<rank_-1,dim,Number>::array_type array_type[dim];

  /**
   * Constructor. Initialize all entries
   * to zero.
   */
  Tensor ();

  /**
   * Copy constructor, where the data is
   * copied from a C-style array.
   */
  Tensor (const array_type &initializer);

  /**
   * Conversion operator from tensor of
   * tensors.
   */
  Tensor (const Tensor<1,dim,Tensor<rank_-1,dim,Number> > &tensor_in);

  /**
   * Conversion operator to tensor of
   * tensors.
   */
  operator Tensor<1,dim,Tensor<rank_-1,dim,Number> > () const;

  /**
   * Read-Write access operator.
   */
  Tensor<rank_-1,dim,Number> &operator [] (const unsigned int i);

  /**
   * Read-only access operator.
   */
  const Tensor<rank_-1,dim,Number> &operator [] (const unsigned int i) const;

  /**
   * Read access using TableIndices <tt>indices</tt>
   */
  Number operator [](const TableIndices<rank_> &indices) const;

  /**
   * Read and write access using TableIndices <tt>indices</tt>
   */
  Number &operator [](const TableIndices<rank_> &indices);

  /**
   *  Assignment operator.
   */
  Tensor &operator = (const Tensor<rank_,dim,Number> &);

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
  Tensor<rank_,dim,Number> &operator = (const Number d);

  /**
   *  Test for equality of two tensors.
   */
  bool operator == (const Tensor<rank_,dim,Number> &) const;

  /**
   *  Test for inequality of two tensors.
   */
  bool operator != (const Tensor<rank_,dim,Number> &) const;

  /**
   *  Add another tensor.
   */
  Tensor<rank_,dim,Number> &operator += (const Tensor<rank_,dim,Number> &);

  /**
   *  Subtract another tensor.
   */
  Tensor<rank_,dim,Number> &operator -= (const Tensor<rank_,dim,Number> &);

  /**
   *  Scale the tensor by <tt>factor</tt>,
   *  i.e. multiply all components by
   *  <tt>factor</tt>.
   */
  Tensor<rank_,dim,Number> &operator *= (const Number factor);

  /**
   *  Scale the vector by
   *  <tt>1/factor</tt>.
   */
  Tensor<rank_,dim,Number> &operator /= (const Number factor);

  /**
   *  Add two tensors. If possible, you
   *  should use <tt>operator +=</tt>
   *  instead since this does not need the
   *  creation of a temporary.
   */
  Tensor<rank_,dim,Number>   operator + (const Tensor<rank_,dim,Number> &) const;

  /**
   *  Subtract two tensors. If possible,
   *  you should use <tt>operator -=</tt>
   *  instead since this does not need the
   *  creation of a temporary.
   */
  Tensor<rank_,dim,Number>   operator - (const Tensor<rank_,dim,Number> &) const;

  /**
   * Unary minus operator. Negate all
   * entries of a tensor.
   */
  Tensor<rank_,dim,Number>   operator - () const;

  /**
   * Return the Frobenius-norm of a tensor,
   * i.e. the square root of the sum of
   * squares of all entries.
   */
  real_type norm () const;

  /**
   * Return the square of the
   * Frobenius-norm of a tensor,
   * i.e. the
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
   * Fill a vector with all tensor elements.
   *
   * This function unrolls all
   * tensor entries into a single,
   * linearly numbered vector. As
   * usual in C++, the rightmost
   * index of the tensor marches fastest.
   */
  template <typename Number2>
  void unroll (Vector<Number2> &result) const;

  /**
   * Returns an unrolled index in
   * the range [0,dim^rank-1] for the element of the tensor indexed by
   * the argument to the function.
   */
  static
  unsigned int
  component_to_unrolled_index(const TableIndices<rank_> &indices);

  /**
   * Opposite of  component_to_unrolled_index: For an index in the
   * range [0,dim^rank-1], return which set of indices it would
   * correspond to.
   */
  static
  TableIndices<rank_> unrolled_to_component_indices(const unsigned int i);



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
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  static std::size_t memory_consumption ();

  /**
   * Exception.
   */
  DeclException1 (ExcInvalidTensorIndex,
                  int,
                  << "Invalid tensor index " << arg1);

  /**
   * Read or write the data of this object to or
   * from a stream for the purpose of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

private:
  /**
   * Array of tensors holding the
   * subelements.
   */
  Tensor<rank_-1,dim,Number> subtensor[dim];

  /**
   * Help function for unroll.
   */
  template <typename Number2>
  void unroll_recursion(Vector<Number2> &result,
                        unsigned int    &start_index) const;

  // make the following class a
  // friend to this class. in principle,
  // it would suffice if otherrank==rank+1,
  // but then the compiler complains
  // that this be an explicit specialization
  // which is not what we want
  //
  // also, it would be sufficient to make
  // the function unroll_loops a friend,
  // but that seems to be impossible as well.
  template <int, int, typename> friend class Tensor;
};


/*--------------------------- Inline functions -----------------------------*/

#ifndef DOXYGEN

template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor ()
{
// default constructor. not specifying an initializer list calls
// the default constructor of the subobjects, which initialize them
// selves. therefore, the tensor is set to zero this way
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] = Tensor<rank_-1,dim,Number>(initializer[i]);
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::Tensor
(const Tensor<1,dim,Tensor<rank_-1,dim,Number> > &tensor_in)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] = tensor_in[i];
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>::operator
Tensor<1,dim,Tensor<rank_-1,dim,Number> > () const
{
  return Tensor<1,dim,Tensor<rank_-1,dim,Number> > (subtensor);
}



template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i)
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  return subtensor[i];
}



template <int rank_, int dim, typename Number>
inline
const typename Tensor<rank_,dim,Number>::value_type &
Tensor<rank_,dim,Number>::operator[] (const unsigned int i) const
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));

  return subtensor[i];
}

template <int rank_, int dim, typename Number>
inline
Number
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices) const
{
  const unsigned int inner_ind = indices[0];
  Assert (inner_ind<dim, ExcIndexRange(inner_ind, 0, dim));

  TableIndices<rank_-1> indices1;
  for (unsigned int i = 0; i < rank_-1; i++)
    indices1[i] = indices[i+1];
  return (subtensor[inner_ind])[indices1];
}

template <int rank_, int dim, typename Number>
inline
Number &
Tensor<rank_,dim,Number>::operator[] (const TableIndices<rank_> &indices)
{
  const unsigned int inner_ind = indices[0];
  Assert (inner_ind<dim, ExcIndexRange(inner_ind, 0, dim));

  TableIndices<rank_-1> indices1;
  for (unsigned int i = 0; i < rank_-1; i++)
    indices1[i] = indices[i+1];
  return (subtensor[inner_ind])[indices1];
}

template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Tensor<rank_,dim,Number> &t)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] = t.subtensor[i];
  return *this;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator = (const Number d)
{
  Assert (d==Number(0), ExcMessage ("Only assignment with zero is allowed"));
  (void) d;

  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] = 0;
  return *this;
}



template <int rank_, int dim, typename Number>
inline
bool
Tensor<rank_,dim,Number>::operator == (const Tensor<rank_,dim,Number> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (subtensor[i] != p.subtensor[i]) return false;
  return true;
}



template <int rank_, int dim, typename Number>
inline
bool
Tensor<rank_,dim,Number>::operator != (const Tensor<rank_,dim,Number> &p) const
{
  return !((*this) == p);
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator += (const Tensor<rank_,dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] += p.subtensor[i];
  return *this;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator -= (const Tensor<rank_,dim,Number> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] -= p.subtensor[i];
  return *this;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator *= (const Number s)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] *= s;
  return *this;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number> &
Tensor<rank_,dim,Number>::operator /= (const Number s)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] /= s;
  return *this;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>
Tensor<rank_,dim,Number>::operator + (const Tensor<rank_,dim,Number> &t) const
{
  Tensor<rank_,dim,Number> tmp(*this);

  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] += t.subtensor[i];

  return tmp;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>
Tensor<rank_,dim,Number>::operator - (const Tensor<rank_,dim,Number> &t) const
{
  Tensor<rank_,dim,Number> tmp(*this);

  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] -= t.subtensor[i];

  return tmp;
}



template <int rank_, int dim, typename Number>
inline
Tensor<rank_,dim,Number>
Tensor<rank_,dim,Number>::operator - () const
{
  Tensor<rank_,dim,Number> tmp;

  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] = -subtensor[i];

  return tmp;
}



template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::real_type
Tensor<rank_,dim,Number>::norm () const
{
  return std::sqrt (norm_square());
}



template <int rank_, int dim, typename Number>
inline
typename Tensor<rank_,dim,Number>::real_type
Tensor<rank_,dim,Number>::norm_square () const
{
  real_type s = 0;
  for (unsigned int i=0; i<dim; ++i)
    s += subtensor[i].norm_square();

  return s;
}



template <int rank_, int dim, typename Number>
template <typename Number2>
inline
void
Tensor<rank_, dim, Number>::unroll (Vector<Number2> &result) const
{
  AssertDimension (result.size(),(Utilities::fixed_power<rank_, unsigned int>(dim)));

  unsigned int index = 0;
  unroll_recursion (result, index);
}



template <int rank_, int dim, typename Number>
template <typename Number2>
inline
void
Tensor<rank_, dim, Number>::unroll_recursion (Vector<Number2> &result,
                                              unsigned int    &index) const
{
  for (unsigned int i=0; i<dim; ++i)
    {
      operator[](i).unroll_recursion(result, index);
    }
}

template <int rank_, int dim, typename Number>
inline
unsigned int
Tensor<rank_, dim, Number>::component_to_unrolled_index(const TableIndices<rank_> &indices)
{
  TableIndices<rank_-1> indices1;
  for (unsigned int i = 0; i < rank_-1; i++)
    indices1[i] = indices[i];

  Assert (indices[rank_-1] < dim,
          ExcIndexRange (indices[rank_-1], 0, dim));
  return ( Tensor<rank_-1,dim,Number>::component_to_unrolled_index(indices1) * dim + indices[rank_-1]);
}

template <int rank_, int dim, typename Number>
inline
TableIndices<rank_>
Tensor<rank_, dim, Number>::unrolled_to_component_indices(const unsigned int i)
{
  Assert (i < n_independent_components,
          ExcIndexRange (i, 0, n_independent_components));

  TableIndices<rank_>   indices;

  unsigned int remainder = i;
  for (int r=rank_-1; r>=0; --r)
    {
      indices[r] = (remainder % dim);
      remainder /= dim;
    }
  Assert (remainder == 0, ExcInternalError());

  return indices;
}



template <int rank_, int dim, typename Number>
inline
void Tensor<rank_,dim,Number>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i].clear();
}



template <int rank_, int dim, typename Number>
inline
std::size_t
Tensor<rank_,dim,Number>::memory_consumption ()
{
  return sizeof(Tensor<rank_,dim,Number>);
}



template <int rank_, int dim, typename Number>
template <class Archive>
inline
void
Tensor<rank_,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &subtensor;
}

#endif // DOXYGEN
/* ----------------- Non-member functions operating on tensors. ------------ */



/**
 * Output operator for tensors. Print the elements consecutively, with
 * a space in between, two spaces between rank 1 subtensors, three
 * between rank 2 and so on.
 *
 * @relates Tensor
 */
template <int rank_, int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out, const Tensor<rank_,dim,Number> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
}

#ifndef DOXYGEN

/**
 * Specialization for 1D.
 */
template <int rank_>
inline
std::ostream &operator << (std::ostream &out, const Tensor<rank_,1> &p)
{
  out << p[0];

  return out;
}

#endif // DOXYGEN


/**
 * Contract a tensor of rank 1 with a tensor of rank 1. The result is
 * <tt>sum_j src1[j] src2[j]</tt>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2000
 */
template <int dim, typename Number>
inline
Number contract (const Tensor<1,dim,Number> &src1,
                 const Tensor<1,dim,Number> &src2)
{
  Number res = 0.;
  for (unsigned int i=0; i<dim; ++i)
    res += src1[i] * src2[i];

  return res;
}


/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Number
operator * (const Tensor<1,dim,Number> &src1,
            const Tensor<1,dim,Number> &src2)
{
  return contract(src1, src2);
}


/**
 * Double contract two tensors of rank 2, thus computing the Frobenius
 * inner product <tt> sum<sub>i,j</sub> src1[i][j]*src2[i][j]</tt>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2000
 */
template <int dim, typename Number>
inline
Number double_contract (const Tensor<2, dim, Number> &src1,
                        const Tensor<2, dim, Number> &src2)
{
  Number res = 0.;
  for (unsigned int i=0; i<dim; ++i)
    res += contract(src1[i],src2[i]);

  return res;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 1. The result is
 * <tt>dest[i] = sum_j src1[i][j] src2[j]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<1,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[i][0] * src2[0];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[i][j] * src2[j];
    }
}


/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
Tensor<1,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<1,dim,Number> &src2)
{
  Tensor<1,dim,Number> dest (false);
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[i][0] * src2[0];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[i][j] * src2[j];
    }
  return dest;
}


/**
 * Contract a tensor of rank 1 with a tensor of rank 2. The result is
 * <tt>dest[i] = sum_j src1[j] src2[j][i]</tt>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void contract (Tensor<1,dim,Number>       &dest,
               const Tensor<1,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[0] * src2[0][i];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[j] * src2[j][i];
    }
}


/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<1,dim,Number>
operator * (const Tensor<1,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<1,dim,Number> dest (false);
  for (unsigned int i=0; i<dim; ++i)
    {
      dest[i] = src1[0] * src2[0][i];
      for (unsigned int j=1; j<dim; ++j)
        dest[i] += src1[j] * src2[j][i];
    }
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The result is
 * <tt>dest[i][k] = sum_j src1[i][j] src2[j][k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
        dest[i][j] = src1[i][0] * src2[0][j];
        for (unsigned int k=1; k<dim; ++k)
          dest[i][j] += src1[i][k] * src2[k][j];
      }
}



/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<2,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dest[i][j] += src1[i][k] * src2[k][j];
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The
 * contraction is performed over index <tt>index1</tt> of the first tensor,
 * and <tt>index2</tt> of the second tensor. Thus, if <tt>index1==2</tt>,
 * <tt>index2==1</tt>, the result is the usual contraction, but if for
 * example <tt>index1==1</tt>, <tt>index2==2</tt>, then the result is
 * <tt>dest[i][k] = sum_j src1[j][i] src2[k][j]</tt>.
 *
 * Note that the number of the index is counted from 1 on, not from
 * zero as usual.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,   const unsigned int index1,
               const Tensor<2,dim,Number> &src2,   const unsigned int index2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[k][i] * src2[k][j];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[k][i] * src2[j][k];
          break;

        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index2)));
        };
      break;
    case 2:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[i][k] * src2[k][j];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                dest[i][j] += src1[i][k] * src2[j][k];
          break;

        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index2)));
        };
      break;

    default:
      Assert (false, (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index1)));
    };
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 1. The
 * contraction is performed over index <tt>index1</tt> of the first
 * tensor.
 *
 * Note that the number of the index is counted from 1 on, not from
 * zero as usual.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<2,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1,   const unsigned int index1,
               const Tensor<1,dim,Number> &src2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[k][i][j] * src2[k];
      break;

    case 2:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[i][k][j] * src2[k];
      break;

    case 3:
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            dest[i][j] += src1[i][j][k] * src2[k];
      break;

    default:
      Assert (false,
              (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index1)));
    };
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The result is
 * <tt>dest[i][j][l] = sum_k src1[i][j][k] src2[k][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1,
               const Tensor<2,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][j][l] * src2[l][k];
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The
 * contraction is performed over index <tt>index1</tt> of the first tensor,
 * and <tt>index2</tt> of the second tensor. Thus, if <tt>index1==3</tt>,
 * <tt>index2==1</tt>, the result is the usual contraction, but if for
 * example <tt>index1==1</tt>, <tt>index2==2</tt>, then the result is
 * <tt>dest[i][j][k] = sum_l src1[l][i][j] src2[k][l]</tt>.
 *
 * Note that the number of the index is counted from 1 on, not from
 * zero as usual.
 *
 * @relates Tensor
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<3,dim,Number> &src1, const unsigned int index1,
               const Tensor<2,dim,Number> &src2, const unsigned int index2)
{
  dest.clear ();

  switch (index1)
    {
    case 1:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[l][i][j] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[l][i][j] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index2)));
        }

      break;
    case 2:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][l][j] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][l][j] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index2)));
        }

      break;
    case 3:
      switch (index2)
        {
        case 1:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][j][l] * src2[l][k];
          break;
        case 2:
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              for (unsigned int k=0; k<dim; ++k)
                for (unsigned int l=0; l<dim; ++l)
                  dest[i][j][k] += src1[i][j][l] * src2[k][l];
          break;
        default:
          Assert (false,
                  (typename Tensor<2,dim,Number>::ExcInvalidTensorIndex (index2)));
        }

      break;
    default:
      Assert (false,
              (typename Tensor<3,dim,Number>::ExcInvalidTensorIndex (index1)));
    }
}


/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<3,dim,Number>
operator * (const Tensor<3,dim,Number> &src1,
            const Tensor<2,dim,Number> &src2)
{
  Tensor<3,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][j][l] * src2[l][k];
  return dest;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 3. The result is
 * <tt>dest[i][j][l] = sum_k src1[i][k] src2[k][j][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
void contract (Tensor<3,dim,Number>       &dest,
               const Tensor<2,dim,Number> &src1,
               const Tensor<3,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][l] * src2[l][j][k];
}


/**
 * Multiplication operator performing a contraction of the last index
 * of the first argument and the first index of the second
 * argument. This function therefore does the same as the
 * corresponding <tt>contract</tt> function, but returns the result as
 * a return value, rather than writing it into the reference given as
 * the first argument to the <tt>contract</tt> function.
 *
 * Note that for the <tt>Tensor</tt> class, the multiplication
 * operator only performs a contraction over a single pair of
 * indices. This is in contrast to the multiplication operator for
 * symmetric tensors, which does the double contraction.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Tensor<3,dim,Number>
operator * (const Tensor<2,dim,Number> &src1,
            const Tensor<3,dim,Number> &src2)
{
  Tensor<3,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j][k] += src1[i][l] * src2[l][j][k];
  return dest;
}


/**
 * Contract a tensor of rank 3 with a tensor of rank 3. The result is
 * <tt>dest[i][j][k][l] = sum_m src1[i][j][m] src2[m][k][l]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Number>
inline
Tensor<4,dim,Number>
operator * (const Tensor<3,dim,Number> &src1,
            const Tensor<3,dim,Number> &src2)
{
  Tensor<4,dim,Number> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          for (unsigned int m=0; m<dim; ++m)
            dest[i][j][k][l] += src1[i][j][m] * src2[m][k][l];
  return dest;
}


/**
 * Contract the last two indices of <tt>src1</tt> with the two indices
 * <tt>src2</tt>, creating a rank-2 tensor. This is the matrix-vector
 * product analog operation between tensors of rank 4 and rank 2.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
void double_contract (Tensor<2,dim,Number>       &dest,
                      const Tensor<4,dim,Number> &src1,
                      const Tensor<2,dim,Number> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          dest[i][j] += src1[i][j][k][l] * src2[k][l];
}


/**
 * Contract three tensors, corresponding to the matrix vector product
 * <i>u<sup>T</sup> A v</i>.
 *
 * @relates Tensor
 * @author Guido Kanschat, 2004
 */
template <int dim, typename Number>
inline
Number contract3 (const Tensor<1,dim,Number> &u,
                  const Tensor<2,dim,Number> &A,
                  const Tensor<1,dim,Number> &v)
{
  Number result = 0.;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      result += u[i] * A[i][j] * v[j];
  return result;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k}
 * a_{i}b_{ijk}c_{jk}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<1,dim,Number> &t1,
           const Tensor<3,dim,Number> &t2,
           const Tensor<2,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        s += t1[i] * t2[i][j][k] * t3[j][k];
  return s;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k}
 * a_{ij}b_{ijk}c_{k}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<2,dim,Number> &t1,
           const Tensor<3,dim,Number> &t2,
           const Tensor<1,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        s += t1[i][j] * t2[i][j][k] * t3[k];
  return s;
}


/**
 * Compute the contraction of three tensors $s=\sum_{i,j,k,l}
 * a_{ij}b_{ijkl}c_{kl}$.
 *
 * @relates Tensor
 * @author Toby D. Young, 2011
 */
template <int dim, typename Number>
inline
Number
contract3 (const Tensor<2,dim,Number> &t1,
           const Tensor<4,dim,Number> &t2,
           const Tensor<2,dim,Number> &t3)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          s += t1[i][j] * t2[i][j][k][l] * t3[k][l];
  return s;
}


/**
 * Form the outer product of two tensors of rank 1 and 1, i.e.
 * <tt>dst[i][j] = src1[i] * src2[j]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<2,dim,Number>       &dst,
                    const Tensor<1,dim,Number> &src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dst[i][j] = src1[i] * src2[j];
}


/**
 * Form the outer product of two tensors of rank 1 and 2, i.e.
 * <tt>dst[i][j][k] = src1[i] * src2[j][k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<3,dim,Number>       &dst,
                    const Tensor<1,dim,Number> &src1,
                    const Tensor<2,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dst[i][j][k] = src1[i] * src2[j][k];
}


/**
 * Form the outer product of two tensors of rank 2 and 1, i.e.
 * <tt>dst[i][j][k] = src1[i][j] * src2[k]</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<3,dim,Number>       &dst,
                    const Tensor<2,dim,Number> &src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        dst[i][j][k] = src1[i][j] * src2[k];
}


/**
 * Form the outer product of two tensors of rank 0 and 1, i.e.
 * <tt>dst[i] = src1 * src2[i]</tt>. Of course, this is only a scaling of
 * <tt>src2</tt>, but we consider this an outer product for completeness of
 * these functions and since this is sometimes needed when writing
 * templates that depend on the rank of a tensor, which may sometimes
 * be zero (i.e. a scalar).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Number                src1,
                    const Tensor<1,dim,Number> &src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1 * src2[i];
}



/**
 * Form the outer product of two tensors of rank 1 and 0, i.e.
 * <tt>dst[i] = src1[i] * src2</tt>. Of course, this is only a scaling of
 * <tt>src1</tt>, but we consider this an outer product for completeness of
 * these functions and since this is sometimes needed when writing
 * templates that depend on the rank of a tensor, which may sometimes
 * be zero (i.e. a scalar).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
void outer_product (Tensor<1,dim,Number>       &dst,
                    const Tensor<1,dim,Number>  src1,
                    const Number         src2)
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1[i] * src2;
}


/**
 * Cross-product in 2d. This is just a rotation by 90 degrees
 * clockwise to compute the outer normal from a tangential
 * vector. This function is defined for all space dimensions to allow
 * for dimension independent programming (e.g. within switches over
 * the space dimenion), but may only be called if the actual dimension
 * of the arguments is two (e.g. from the <tt>dim==2</tt> case in the
 * switch).
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src)
{
  Assert (dim==2, ExcInternalError());

  dst[0] = src[1];
  dst[1] = -src[0];
}


/**
 * Cross-product of 2 vectors in 3d. This function is defined for all
 * space dimensions to allow for dimension independent programming
 * (e.g. within switches over the space dimenion), but may only be
 * called if the actual dimension of the arguments is three (e.g. from
 * the <tt>dim==3</tt> case in the switch).
 *
 * @relates Tensor
 * @author Guido Kanschat, 2001
 */
template <int dim, typename Number>
inline
void
cross_product (Tensor<1,dim,Number>       &dst,
               const Tensor<1,dim,Number> &src1,
               const Tensor<1,dim,Number> &src2)
{
  Assert (dim==3, ExcInternalError());

  dst[0] = src1[1]*src2[2] - src1[2]*src2[1];
  dst[1] = src1[2]*src2[0] - src1[0]*src2[2];
  dst[2] = src1[0]*src2[1] - src1[1]*src2[0];
}


/**
 * Compute the scalar product $a:b=\sum_{i,j} a_{ij}b_{ij}$ between two
 * tensors $a,b$ of rank 2. We don't use <code>operator*</code> for this
 * operation since the product between two tensors is usually assumed to be
 * the contraction over the last index of the first tensor and the first index
 * of the second tensor, for example $(a\cdot b)_{ij}=\sum_k a_{ik}b_{kj}$.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2008
 */
template <int dim, typename Number>
inline
Number
scalar_product (const Tensor<2,dim,Number> &t1,
                const Tensor<2,dim,Number> &t2)
{
  Number s = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      s += t1[i][j] * t2[i][j];
  return s;
}


/**
 * Compute the determinant of a tensor of arbitrary rank and dimension
 * one. Since this is a number, the return value is, of course, the
 * number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <int rank, typename Number>
inline
Number determinant (const Tensor<rank,1,Number> &t)
{
  return determinant(t[0]);
}



/**
 * Compute the determinant of a tensor of rank one and dimension
 * one. Since this is a number, the return value is, of course, the
 * number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<1,1,Number> &t)
{
  return t[0];
}


/**
 * Compute the determinant of a tensor of rank two and dimension
 * one. Since this is a number, the return value is, of course, the
 * number itself.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,1,Number> &t)
{
  return t[0][0];
}



/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==2</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,2,Number> &t)
{
  return ((t[0][0] * t[1][1]) - (t[1][0] * t[0][1]));
}


/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 1998
 */
template <typename Number>
inline
Number determinant (const Tensor<2,3,Number> &t)
{
  // use exactly the same expression with the
  // same order of operations as for the inverse
  // to let the compiler use common
  // subexpression elimination when using
  // determinant and inverse in nearby code
  const Number t4 = t[0][0]*t[1][1],
               t6 = t[0][0]*t[1][2],
               t8 = t[0][1]*t[1][0],
               t00 = t[0][2]*t[1][0],
               t01 = t[0][1]*t[2][0],
               t04 = t[0][2]*t[2][0],
               det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                      t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
  return det;
}


/**
 * Compute the determinant of a tensor or rank 2, here for all dimensions for
 * which no explicit specialization is available above.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2009
 */
template <int dim, typename Number>
inline
Number determinant (const Tensor<2,dim,Number> &t)
{
  // compute the determinant using the
  // Laplace expansion of the
  // determinant. this may not be the most
  // efficient algorithm, but it does for
  // small n.
  //
  // for some algorithmic simplicity, we use
  // the expansion along the last row
  Number det = 0;

  for (unsigned int k=0; k<dim; ++k)
    {
      Tensor<2,dim-1,Number> minor;
      for (unsigned int i=0; i<dim-1; ++i)
        for (unsigned int j=0; j<dim-1; ++j)
          minor[i][j] = t[i][j<k ? j : j+1];

      const Number cofactor = std::pow (-1., static_cast<Number>(k+1)) *
                              determinant (minor);

      det += t[dim-1][k] * cofactor;
    }

  return std::pow (-1., static_cast<Number>(dim)) * det;
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of
 * its diagonal entries.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2001
 */
template <int dim, typename Number>
Number trace (const Tensor<2,dim,Number> &d)
{
  Number t=d[0][0];
  for (unsigned int i=1; i<dim; ++i)
    t += d[i][i];
  return t;
}



/**
 * Compute and return the inverse of the given tensor. Since the
 * compiler can perform the return value optimization, and since the
 * size of the return object is known, it is acceptable to return the
 * result by value, rather than by reference as a parameter.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
invert (const Tensor<2,dim,Number> &t)
{
  Number return_tensor [dim][dim];
  switch (dim)
    {
    case 1:
      return_tensor[0][0] = 1.0/t[0][0];
      break;

    case 2:
      // this is Maple output,
      // thus a bit unstructured
    {
      const Number det = t[0][0]*t[1][1]-t[1][0]*t[0][1];
      const Number t4 = 1.0/det;
      return_tensor[0][0] = t[1][1]*t4;
      return_tensor[0][1] = -t[0][1]*t4;
      return_tensor[1][0] = -t[1][0]*t4;
      return_tensor[1][1] = t[0][0]*t4;
      break;
    }

    case 3:
    {
      const Number t4 = t[0][0]*t[1][1],
                   t6 = t[0][0]*t[1][2],
                   t8 = t[0][1]*t[1][0],
                   t00 = t[0][2]*t[1][0],
                   t01 = t[0][1]*t[2][0],
                   t04 = t[0][2]*t[2][0],
                   det = (t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
                          t00*t[2][1]+t01*t[1][2]-t04*t[1][1]),
                         t07 = 1.0/det;
      return_tensor[0][0] = (t[1][1]*t[2][2]-t[1][2]*t[2][1])*t07;
      return_tensor[0][1] = (t[0][2]*t[2][1]-t[0][1]*t[2][2])*t07;
      return_tensor[0][2] = (t[0][1]*t[1][2]-t[0][2]*t[1][1])*t07;
      return_tensor[1][0] = (t[1][2]*t[2][0]-t[1][0]*t[2][2])*t07;
      return_tensor[1][1] = (t[0][0]*t[2][2]-t04)*t07;
      return_tensor[1][2] = (t00-t6)*t07;
      return_tensor[2][0] = (t[1][0]*t[2][1]-t[1][1]*t[2][0])*t07;
      return_tensor[2][1] = (t01-t[0][0]*t[2][1])*t07;
      return_tensor[2][2] = (t4-t8)*t07;

      break;
    }

    // if desired, take over the
    // inversion of a 4x4 tensor
    // from the FullMatrix
    default:
      AssertThrow (false, ExcNotImplemented());
    }
  return Tensor<2,dim,Number>(return_tensor);
}



/**
 * Return the transpose of the given tensor. Since the compiler can
 * perform the return value optimization, and since the size of the
 * return object is known, it is acceptable to return the result by
 * value, rather than by reference as a parameter. Note that there are
 * specializations of this function for <tt>dim==1,2,3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <int dim, typename Number>
inline
Tensor<2,dim,Number>
transpose (const Tensor<2,dim,Number> &t)
{
  Number tt[dim][dim];
  for (unsigned int i=0; i<dim; ++i)
    {
      tt[i][i] = t[i][i];
      for (unsigned int j=i+1; j<dim; ++j)
        {
          tt[i][j] = t[j][i];
          tt[j][i] = t[i][j];
        };
    }
  return Tensor<2,dim,Number>(tt);
}

#ifndef DOXYGEN

/**
 * Return the transpose of the given tensor. This is the
 * specialization of the general template for <tt>dim==1</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,1,Number>
transpose (const Tensor<2,1,Number> &t)
{
  return t;
}




/**
 * Return the transpose of the given tensor. This is the
 * specialization of the general template for <tt>dim==2</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,2,Number>
transpose (const Tensor<2,2,Number> &t)
{
  const Number x[2][2] = {{t[0][0], t[1][0]}, {t[0][1], t[1][1]}};
  return Tensor<2,2,Number>(x);
}




/**
 * Return the transpose of the given tensor. This is the
 * specialization of the general template for <tt>dim==3</tt>.
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2002
 */
template <typename Number>
inline
Tensor<2,3,Number>
transpose (const Tensor<2,3,Number> &t)
{
  const Number x[3][3] = {{t[0][0], t[1][0], t[2][0]},
    {t[0][1], t[1][1], t[2][1]},
    {t[0][2], t[1][2], t[2][2]}
  };
  return Tensor<2,3,Number>(x);
}

#endif // DOXYGEN


/**
 * Return the $l_1$ norm of the given rank-2 tensor, where
 * $||t||_1 = \max_j \sum_i |t_{ij}|$ (maximum of
 * the sums over columns).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
double
l1_norm (const Tensor<2,dim,Number> &t)
{
  double max = 0;
  for (unsigned int j=0; j<dim; ++j)
    {
      double sum = 0;
      for (unsigned int i=0; i<dim; ++i)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}



/**
 * Return the $l_\infty$ norm of the given rank-2 tensor, where
 * $||t||_\infty = \max_i \sum_j |t_{ij}|$ (maximum of
 * the sums over rows).
 *
 * @relates Tensor
 * @author Wolfgang Bangerth, 2012
 */
template <int dim, typename Number>
inline
double
linfty_norm (const Tensor<2,dim,Number> &t)
{
  double max = 0;
  for (unsigned int i=0; i<dim; ++i)
    {
      double sum = 0;
      for (unsigned int j=0; j<dim; ++j)
        sum += std::fabs(t[i][j]);

      if (sum > max)
        max = sum;
    }

  return max;
}



/**
 * Multiplication of a tensor of general rank with a scalar Number
 * from the right.
 *
 * @relates Tensor
 */
template <int rank, int dim, typename Number>
inline
Tensor<rank,dim,Number>
operator * (const Tensor<rank,dim,Number> &t,
            const Number                   factor)
{
  Tensor<rank,dim,Number> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a tensor of general rank with a scalar Number
 * from the left.
 *
 * @relates Tensor
 */
template <int rank, int dim, typename Number>
inline
Tensor<rank,dim,Number>
operator * (const Number                   factor,
            const Tensor<rank,dim,Number> &t)
{
  Tensor<rank,dim,Number> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Division of a tensor of general rank by a scalar Number.
 *
 * @relates Tensor
 */
template <int rank, int dim, typename Number>
inline
Tensor<rank,dim,Number>
operator / (const Tensor<rank,dim,Number> &t,
            const Number                   factor)
{
  Tensor<rank,dim,Number> tt = t;
  tt /= factor;
  return tt;
}




/**
 * Multiplication of a tensor of general rank with a scalar double
 * from the right.
 *
 * @relates Tensor
 */
template <int rank, int dim>
inline
Tensor<rank,dim>
operator * (const Tensor<rank,dim> &t,
            const double            factor)
{
  Tensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a tensor of general rank with a scalar double
 * from the left.
 *
 * @relates Tensor
 */
template <int rank, int dim>
inline
Tensor<rank,dim>
operator * (const double            factor,
            const Tensor<rank,dim> &t)
{
  Tensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Division of a tensor of general rank by a scalar double.
 *
 * @relates Tensor
 */
template <int rank, int dim>
inline
Tensor<rank,dim>
operator / (const Tensor<rank,dim> &t,
            const double            factor)
{
  Tensor<rank,dim> tt = t;
  tt /= factor;
  return tt;
}



/**
 * Multiplication of a tensor of general rank by a scalar
 * complex<double> from the left.
 *
 * @relates Tensor
 */
template <int rank, int dim>
inline
Tensor<rank,dim,std::complex<double> >
operator * (const std::complex<double>  factor,
            const Tensor<rank,dim>     &t)
{
  Tensor<rank,dim,std::complex<double> > tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = factor * t[d];
  return tt;
}



/**
 * Multiplication of a tensor of general rank by a scalar
 * complex<double> from the right.
 *
 * @relates Tensor
 */
template <int rank, int dim>
inline
Tensor<rank,dim,std::complex<double> >
operator * (const Tensor<rank,dim>     &t,
            const std::complex<double>  factor)
{
  Tensor<rank,dim,std::complex<double> > tt;
  for (unsigned int d=0; d<dim; ++d)
    tt[d] = t[d] * factor;
  return tt;
}



DEAL_II_NAMESPACE_CLOSE

#endif
