//----------------------------  tensor.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor.h  ---------------------------
#ifndef __deal2__tensor_h
#define __deal2__tensor_h


#include <base/tensor_base.h>

template <int rank_, int dim> class Tensor;

/**
 * Provide a general tensor class with an arbitrary rank, i.e. with
 * an arbitrary number of indices. The Tensor class provides an
 * indexing operator and a bit of infrastructure, but most
 * functionality is recursively handed down to tensors of rank 1 or
 * put into external templated functions, e.g. the @p{contract} family.
 *
 * Using this tensor class for objects of rank 2 has advantages over
 * matrices in many cases since the dimension is known to the compiler
 * as well as the location of the data. It is therefore possible to
 * produce far more efficient code than for matrices with
 * runtime-dependent dimension.
 */
template <int rank_, int dim>
class Tensor
{
  public:
				     /**
				      * Provide a way to get the
				      * dimension of an object without
				      * explicit knowledge of it's
				      * data type. Implementation is
				      * this way instead of providing
				      * a function @p{dimension()}
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
				      * Declare an array type which can
				      * be used to initialize statically
				      * an object of this type.
				      */
    typedef typename Tensor<rank_-1,dim>::array_type array_type[dim];

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
				      * Read-Write access operator.
				      */
    Tensor<rank_-1,dim> &operator [] (const unsigned int i);

				     /**
				      * Read-only access operator.
				      */
    const Tensor<rank_-1,dim> &operator [] (const unsigned int i) const;

				     /**
				      *  Assignment operator.
				      */
    Tensor & operator = (const Tensor<rank_,dim> &);

				     /**
				      *  Test for equality of two tensors.
				      */
    bool operator == (const Tensor<rank_,dim> &) const;

    				     /**
				      *  Test for inequality of two tensors.
				      */
    bool operator != (const Tensor<rank_,dim> &) const;

				     /**
				      *  Add another vector, i.e. move this tensor by
				      *  the given offset.
				      */
    Tensor<rank_,dim> & operator += (const Tensor<rank_,dim> &);
    
				     /**
				      *  Subtract another vector.
				      */
    Tensor<rank_,dim> & operator -= (const Tensor<rank_,dim> &);

				     /**
				      *  Scale the tensor by @p{factor}, i.e. multiply
				      *  all coordinates by @p{factor}.
				      */
    Tensor<rank_,dim> & operator *= (const double &factor);

				     /**
				      *  Scale the vector by @p{1/factor}.
				      */
    Tensor<rank_,dim> & operator /= (const double &factor);

				     /**
				      *  Add two tensors. If possible, use
				      *  @p{operator +=} instead since this does not
				      *  need to copy a point at least once.
				      */
    Tensor<rank_,dim>   operator + (const Tensor<rank_,dim> &) const;

				     /**
				      *  Subtract two tensors. If possible, use
				      *  @p{operator +=} instead since this does not
				      *  need to copy a point at least once.
				      */
    Tensor<rank_,dim>   operator - (const Tensor<rank_,dim> &) const;

				     /**
				      * Invert all entries of a tensor.
				      */
    Tensor<rank_,dim>   operator - () const;
    
				     /**
				      * Fill a vector with all tensor elements.
				      *
				      * This function unrolls all
				      * tensor entries into a single,
				      * linearly numbered vector. As
				      * usual in C++, the rightmost
				      * index of the tensor marches fastest.
				      */
    void unroll(Vector<double> & result) const;


				     /**
				      * Reset all values to zero.
				      */
    void clear ();

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object.
				      */
    static unsigned int memory_consumption ();

// 				     /**
// 				      * Exception
// 				      */
//     DeclException2(ExcWrongVectorSize, unsigned, int, << "Tensor has " << arg1
// 		   << " entries, but vector has size " << arg2);
    
  private:
				     /**
				      * Array of tensors holding the
				      * subelements.
				      */
    Tensor<rank_-1,dim> subtensor[dim];

				     /**
				      * Help function for unroll.
				      */
    void unroll_recursion(Vector<double> &result,
			  unsigned int   &start_index) const;

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
    template<int otherrank, int otherdim> friend class Tensor;
};


/*--------------------------- Inline functions -----------------------------*/



template <int rank_, int dim>
inline
Tensor<rank_,dim>::Tensor ()
{
// default constructor. not specifying an initializer list calls
// the default constructor of the subobjects, which initialize them
// selves. therefore, the tensor is set to zero this way
};


template <int rank_, int dim>
inline
Tensor<rank_,dim>::Tensor (const array_type &initializer)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] =  Tensor<rank_-1,dim>(initializer[i]);
};


template <int rank_, int dim>
inline
Tensor<rank_-1,dim> &
Tensor<rank_,dim>::operator[] (const unsigned int i)
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  
  return subtensor[i];
};


template <int rank_, int dim>
inline
const Tensor<rank_-1,dim> &
Tensor<rank_,dim>::operator[] (const unsigned int i) const
{
  Assert (i<dim, ExcIndexRange(i, 0, dim));
  
  return subtensor[i];
};


template <int rank_, int dim>
inline
Tensor<rank_,dim> &
Tensor<rank_,dim>::operator = (const Tensor<rank_,dim> &t)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] = t.subtensor[i];
  return *this;
};


template <int rank_, int dim>
inline
bool
Tensor<rank_,dim>::operator == (const Tensor<rank_,dim> &p) const
{
  for (unsigned int i=0; i<dim; ++i)
    if (subtensor[i] != p.subtensor[i]) return false;
  return true;
};


template <int rank_, int dim>
inline
bool
Tensor<rank_,dim>::operator != (const Tensor<rank_,dim> &p) const
{
  return !((*this) == p);
};


template <int rank_, int dim>
inline
Tensor<rank_,dim> &
Tensor<rank_,dim>::operator += (const Tensor<rank_,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] += p.subtensor[i];
  return *this;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim> &
Tensor<rank_,dim>::operator -= (const Tensor<rank_,dim> &p)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] -= p.subtensor[i];
  return *this;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim> &
Tensor<rank_,dim>::operator *= (const double &s)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] *= s;
  return *this;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim> &
Tensor<rank_,dim>::operator /= (const double &s)
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i] /= s;
  return *this;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim>
Tensor<rank_,dim>::operator + (const Tensor<rank_,dim> &t) const
{
  Tensor<rank_,dim> tmp(*this);
  
  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] += t.subtensor[i];

  return tmp;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim>
Tensor<rank_,dim>::operator - (const Tensor<rank_,dim> &t) const
{
  Tensor<rank_,dim> tmp(*this);
  
  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] -= t.subtensor[i];

  return tmp;
};


template <int rank_, int dim>
inline
Tensor<rank_,dim>
Tensor<rank_,dim>::operator - () const
{
  Tensor<rank_,dim> tmp;
  
  for (unsigned int i=0; i<dim; ++i)
    tmp.subtensor[i] = -subtensor[i];

  return tmp;
};


template <int rank_, int dim>
inline
void Tensor<rank_,dim>::clear ()
{
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i].clear();
};



template <int rank_, int dim>
inline
unsigned int
Tensor<rank_,dim>::memory_consumption ()
{
  return sizeof(Tensor<rank_,dim>);
};


/* ----------------- Non-member functions operating on tensors. ------------ */


/*  Exception class. This is certainly not the best possible place for its
    declaration, but at present, local classes to any of Tensor<> can't
    be properly accessed (haven't investigated why). If anyone has a better
    idea, realize it!
*/
DeclException1 (ExcInvalidTensorIndex,
		int,
		<< "Invalid tensor index " << arg1);


/**
 * Output operator for tensors. Print the elements consecutively, with
 * a space in between, two spaces between rank 1 subtensors, three
 * between rank 2 and so on.
 */
template <int rank_, int dim>
inline
std::ostream & operator << (std::ostream &out, const Tensor<rank_,dim> &p)
{
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  return out;
};


/**
 * Specialization for 1D.
 */
template <int rank_>
inline
std::ostream & operator << (std::ostream &out, const Tensor<rank_,1> &p)
{
  out << p[0];

  return out;
};



/**
 * Contract a tensor of rank 1 with a tensor of rank 1. The result is
 * @p{sum_j src1[j] src2[j]}.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
inline
double contract (const Tensor<1,dim> &src1,
		 const Tensor<1,dim> &src2)
{
  double res = 0.;
  for (unsigned int i=0; i<dim; ++i)
    res += src1[i] * src2[i];

  return res;
}


/**
 * Contract a tensor of rank 2 with a tensor of rank 1. The result is
 * @p{dest[i] = sum_j src1[i][j] src2[j]}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<1,dim>       &dest,
	       const Tensor<2,dim> &src1,
	       const Tensor<1,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dest[i] += src1[i][j] * src2[j];
};



/**
 * Contract a tensor of rank 1 with a tensor of rank 2. The result is
 * @p{dest[i] = sum_j src1[j] src2[j][i]}.
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
inline
void contract (Tensor<1,dim>       &dest,
	       const Tensor<1,dim> &src1,
	       const Tensor<2,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dest[i] += src1[j] * src2[j][i];
};



/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The result is
 * @p{dest[i][k] = sum_j src1[i][j] src2[j][k]}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<2,dim>       &dest,
	       const Tensor<2,dim> &src1,
	       const Tensor<2,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	dest[i][j] += src1[i][k] * src2[k][j];
};



/**
 * Contract a tensor of rank 2 with a tensor of rank 2. The
 * contraction is performed over index @p{index1} of the first tensor,
 * and @p{index2} of the second tensor. Thus, if @p{index1==2},
 * @p{index2==1}, the result is the usual contraction, but if for
 * example @p{index1==1}, @p{index2==2}, then the result is
 * @p{dest[i][k] = sum_j src1[j][i] src2[k][j]}.
 *
 * Note that the number of the index is counted from 1 on, not from
 * zero as usual.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<2,dim>       &dest,
	       const Tensor<2,dim> &src1,   const unsigned int index1,
	       const Tensor<2,dim> &src2,   const unsigned int index2)
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
		      Assert (false, ExcInvalidTensorIndex (index2));
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
		      Assert (false, ExcInvalidTensorIndex (index2));
	      };
	    break;

      default:
	    Assert (false, ExcInvalidTensorIndex (index1));
    };
};



/**
 * Contract a tensor of rank 3 with a tensor of rank 1. The
 * contraction is performed over index @p{index1} of the first
 * tensor.
 *
 * Note that the number of the index is counted from 1 on, not from
 * zero as usual.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<2,dim>       &dest,
	       const Tensor<3,dim> &src1,   const unsigned int index1,
	       const Tensor<1,dim> &src2)
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
	    Assert (false, ExcInvalidTensorIndex (index1));
    };
};



/**
 * Contract a tensor of rank 3 with a tensor of rank 2. The result is
 * @p{dest[i][j][l] = sum_k src1[i][j][k] src2[k][l]}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<3,dim>       &dest,
	       const Tensor<3,dim> &src1,
	       const Tensor<2,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  dest[i][j][k] += src1[i][j][l] * src2[l][k];
};



/**
 * Contract a tensor of rank 2 with a tensor of rank 3. The result is
 * @p{dest[i][j][l] = sum_k src1[i][k] src2[k][j][l]}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<3,dim>       &dest,
	       const Tensor<2,dim> &src1,
	       const Tensor<3,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  dest[i][j][k] += src1[i][l] * src2[l][j][k];
};


/**
 * Contract a tensor of rank 3 with a tensor of rank 3. The result is
 * @p{dest[i][j][k][l] = sum_m src1[i][j][m] src2[m][k][l]}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
inline
void contract (Tensor<4,dim>       &dest,
	       const Tensor<3,dim> &src1,
	       const Tensor<3,dim> &src2)
{
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  for (unsigned int m=0; m<dim; ++m)
	    dest[i][j][k][l] += src1[i][j][m] * src2[m][k][l];
};



/**
 * Form the outer product of two tensors of rank 1 and 1, i.e.
 * @p{dst[i][j] = src1[i] * src2[j]}.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
void outer_product (Tensor<2,dim>       &dst,
		    const Tensor<1,dim> &src1,
		    const Tensor<1,dim> &src2) 
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dst[i][j] = src1[i] * src2[j];
};



/**
 * Form the outer product of two tensors of rank 1 and 2, i.e.
 * @p{dst[i][j][k] = src1[i] * src2[j][k]}.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
void outer_product (Tensor<3,dim>       &dst,
		    const Tensor<1,dim> &src1,
		    const Tensor<2,dim> &src2) 
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	dst[i][j][k] = src1[i] * src2[j][k];
};



/**
 * Form the outer product of two tensors of rank 2 and 1, i.e.
 * @p{dst[i][j][k] = src1[i][j] * src2[k]}.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
void outer_product (Tensor<3,dim>       &dst,
		    const Tensor<2,dim> &src1,
		    const Tensor<1,dim> &src2) 
{
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	dst[i][j][k] = src1[i][j] * src2[k];
};



/**
 * Form the outer product of two tensors of rank 0 and 1, i.e.
 * @p{dst[i] = src1 * src2[i]}. Of course, this is only a scaling of
 * @p{src2}, but we consider this an outer product for completeness of
 * these functions and since this is sometimes needed when writing
 * templates that depend on the rank of a tensor, which may sometimes
 * be zero (i.e. a scalar).
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
void outer_product (Tensor<1,dim>       &dst,
		    const double         src1,
		    const Tensor<1,dim> &src2) 
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1 * src2[i];
};



/**
 * Form the outer product of two tensors of rank 1 and 0, i.e.
 * @p{dst[i] = src1[i] * src2}. Of course, this is only a scaling of
 * @p{src1}, but we consider this an outer product for completeness of
 * these functions and since this is sometimes needed when writing
 * templates that depend on the rank of a tensor, which may sometimes
 * be zero (i.e. a scalar).
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
void outer_product (Tensor<1,dim>       &dst,
		    const Tensor<1,dim>  src1,
		    const double         src2) 
{
  for (unsigned int i=0; i<dim; ++i)
    dst[i] = src1[i] * src2;
};



/**
 * Cross-product in 2d. This is just a rotation by 90 degrees
 * clockwise to compute the outer normal from a tangential
 * vector. This function is defined for all space dimensions to allow
 * for dimension independent programming (e.g. within switches over
 * the space dimenion), but may only be called if the actual dimension
 * of the arguments is two (e.g. from the @p{dim==2} case in the
 * switch).
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
inline
void
cross_product (Tensor<1,dim>       &dst,
	       const Tensor<1,dim> &src)
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
 * the @p{dim==3} case in the switch).
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
inline
void
cross_product (Tensor<1,dim>       &dst,
	       const Tensor<1,dim> &src1,
	       const Tensor<1,dim> &src2)
{
  Assert (dim==3, ExcInternalError());
  
  dst[0] = src1[1]*src2[2] - src1[2]*src2[1];
  dst[1] = src1[2]*src2[0] - src1[0]*src2[2];
  dst[2] = src1[0]*src2[1] - src1[1]*src2[0];
}



/**
 * Compute the determinant of a tensor of arbitrary rank and dimension
 * one. Since this is a number, the return value is, of course, the
 * number itself.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int rank>
inline
double determinant (const Tensor<rank,1> &t)
{
				   // determinant of tensors of
				   // dimension one and arbitrary rank
				   // can be computed by recursion. we
				   // need therefore not try to access
				   // the number itself, which is
				   // difficult since it needs @p{rank}
				   // indirections, which is not
				   // computable in the general
				   // template
  return determinant(t[0]);
};



/**
 * Compute the determinant of a tensor of rank one and dimension
 * one. Since this is a number, the return value is, of course, the
 * number itself.
 *
 * @author Wolfgang Bangerth, 1998
 */
inline
double determinant (const Tensor<1,1> &t)
{
  return t[0];
};



/**
 * Compute the determinant of a tensor or rank 2, here for @p{dim==2}.
 *
 * @author Wolfgang Bangerth, 1998
 */
inline
double determinant (const Tensor<2,2> &t)
{
  return ((t[0][0] * t[1][1]) -
	  (t[1][0] * t[0][1]));
};




/**
 * Compute the determinant of a tensor or rank 2, here for @p{dim==3}.
 *
 * @author Wolfgang Bangerth, 1998
 */
inline
double determinant (const Tensor<2,3> &t)
{
				   // get this using Maple:
				   // with(linalg);
				   // a := matrix(3,3);
				   // x := det(a);
				   // readlib(C);
				   // C(x, optimized);
  return ( t[0][0]*t[1][1]*t[2][2]
	   -t[0][0]*t[1][2]*t[2][1]
	   -t[1][0]*t[0][1]*t[2][2]
	   +t[1][0]*t[0][2]*t[2][1]
	   +t[2][0]*t[0][1]*t[1][2]
	   -t[2][0]*t[0][2]*t[1][1] );
};



/**
 * Compute and return the inverse of the given tensor. Since the
 * compiler can perform the return value optimization, and since the
 * size of the return object is known, it is acceptable to return the
 * result by value, rather than by reference as a parameter.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
inline
Tensor<2,dim>
invert (const Tensor<2,dim> &t)
{
  Tensor<2,dim> return_tensor;
  switch (dim) 
    {
      case 1:
	    return_tensor[0][0] = 1.0/t[0][0];
	    return return_tensor;
      case 2:
					     // this is Maple output,
					     // thus a bit unstructured
      {
	const double t4 = 1.0/(t[0][0]*t[1][1]-t[0][1]*t[1][0]);
	return_tensor[0][0] = t[1][1]*t4;
	return_tensor[0][1] = -t[0][1]*t4;
	return_tensor[1][0] = -t[1][0]*t4;
	return_tensor[1][1] = t[0][0]*t4;
	return return_tensor;
      };
      
      case 3:
      {
	const double t4 = t[0][0]*t[1][1],
		     t6 = t[0][0]*t[1][2],
		     t8 = t[0][1]*t[1][0],
		    t00 = t[0][2]*t[1][0],
		    t01 = t[0][1]*t[2][0],
		    t04 = t[0][2]*t[2][0],
		    t07 = 1.0/(t4*t[2][2]-t6*t[2][1]-t8*t[2][2]+
			       t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
	return_tensor[0][0] = (t[1][1]*t[2][2]-t[1][2]*t[2][1])*t07;
	return_tensor[0][1] = -(t[0][1]*t[2][2]-t[0][2]*t[2][1])*t07;
	return_tensor[0][2] = -(-t[0][1]*t[1][2]+t[0][2]*t[1][1])*t07;
	return_tensor[1][0] = -(t[1][0]*t[2][2]-t[1][2]*t[2][0])*t07;
	return_tensor[1][1] = (t[0][0]*t[2][2]-t04)*t07;
	return_tensor[1][2] = -(t6-t00)*t07;
	return_tensor[2][0] = -(-t[1][0]*t[2][1]+t[1][1]*t[2][0])*t07;
	return_tensor[2][1] = -(t[0][0]*t[2][1]-t01)*t07;
	return_tensor[2][2] = (t4-t8)*t07;
	return return_tensor;
      };

					// if desired, take over the
					// inversion of a 4x4 tensor
					// from the FullMatrix
       
      default:
	    AssertThrow (false, ExcNotImplemented());
    };    
  return return_tensor;
};


#endif
